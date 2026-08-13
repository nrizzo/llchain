#ifndef VERSION
#define VERSION "GitHub"
#endif

module;
#include <chrono> // fix for GCC 15.2.1

import <iostream>;
import <string>;
import <vector>;
import <fstream>;
import <cassert>;
import utils;
import algo;
import kseq;
import mummer_essaMEM_wrapper;
import htslib_wrapper;
import chainx;
import cli11;

using std::cout, std::cerr, std::endl;
using std::string, std::to_string;
using std::vector;
using std::ifstream, std::ofstream;
using namespace llchain;
using utils::anchor_t, utils::random_anchors, utils::plot_gap_gap_lower_diag, utils::plot_anchors, utils::Image, utils::place_dummy_anchors, utils::weak_sort_anchors, utils::chainx_sort_anchors, utils::merge_perfect_chains, utils::read_mummer_anchors_single, utils::cumulative_length;
typedef std::size_t size_type;

enum anchor_type_e { MUM, MEM };

struct seeding_params { anchor_type_e anchor_type; int anchor_length; };
struct chaining_params {
	anchor_type_e anchortype;
	int anchor_length;
	utils::chaining_mode mode;
	bool chainx;
	bool chainx_optimal;
	chainx::mode chainx_mode;
	bool chainx_original_magic_numbers;
	bool chainx_optimal_ensure_pred;
	bool store_sam_sequence;
};
struct io_streams {
	ifstream custom_anchors_fs;
	ofstream sam_out;
	ofstream paf_out;
	ofstream mummer_out;
	ofstream phylip_out;
};

/*
 * seed and compute the anchored edit distance
 */
utils::anchor_index_t seed_and_chain(
	const seeding_params &seed_p,
	const chaining_params &chain_p,
	io_streams &ioss,
	const mummer_essaMEM_wrapper::sparseSA &index,
	const string &text,
	const string &text_id,
	string &query,
	const string &query_id,
	const bool backtrack, // find the optimal chain, and log the time spent
	const bool is_query_rev_comp // reverse complemented query
);

int main(int argc, char **argv)
{
	seeding_params seed_p;
	chaining_params chain_p;
	io_streams ioss;

	// ### CLI setup ###
	CLI::App app("Log-linear chaining to compute the anchored edit distance", "llchain");
	app.usage("Usage: llchain [-t ref.fa] [-q queries.fa] [-a {MUM,MEM}] [-m {global,semiglobal}] [OPTIONS]");
	app.set_version_flag("--version", "llchain version " + string(VERSION));
	string text_path {""}, query_path {""};
	auto t_opt = app.add_option("-t,--text", text_path, "One or more DNA reference sequences in FASTA format (.fa,.fa.gz)")
		->check(CLI::ExistingFile);
	auto q_opt = app.add_option("-q,--query,--queries", query_path, "DNA queries in FASTA format (.fa,.fa.gz)")
		->check(CLI::ExistingFile);

	string anchor_type_s {}; // to be replaced by seed_p.anchor_type
	auto anchor_type_opt = app.add_option("-a,--anchor-type", anchor_type_s, "MUM or MEM (computed via MUMmer/essaMEM)")
		->default_val("MUM");

	auto anchor_length_opt = app.add_option("-l,--length", seed_p.anchor_length, "Minimum anchor length")
		->check(CLI::Range(2, std::numeric_limits<int>::max()))
		->default_val(20);

	string custom_anchors_path {""};
	app.add_option("--custom-anchors", custom_anchors_path, "Do not index/query but read the anchors from this file (NB it should respect the same order as query file)")
		->check(CLI::ExistingFile)
		->excludes(anchor_type_opt)
		->excludes(anchor_length_opt);

	string mode_s; // to be replaced by chain_p.mode and chain_p.chainx_mode
	app.add_option("-m,--mode", mode_s, "Chaining mode (global or semiglobal)")
		->default_val("global");

	bool reverse_complement {false};
	app.add_flag("-r,--reverse-complement", reverse_complement, "Compute the anchored edit distance for the reversed-and-complemented queries as well");

	bool all_to_all {false};
	auto all_to_all_option = app.add_flag("--all-to-all", all_to_all, "Pairwise comparisons (of the text sequences, see also --phylip)")
		->excludes(q_opt)
		->needs(t_opt);

	auto chainx_option = app.add_flag("--chainx", chain_p.chainx, "Chain with at-cg/ChainX algorithm (variable B >= 100, alpha = 4)")
		->group("")
		->default_val(false);

	auto chainx_optimal_option = app.add_flag("--chainx-opt", chain_p.chainx_optimal, "Chain with algbio/ChainX *optimal* algorithm (variable B >= 100, alpha = 4)")
		->group("")
		->excludes(chainx_option)
		->default_val(false);

	chain_p.chainx_original_magic_numbers = false;
	app.add_flag("--chainx-original-magic-numbers", chain_p.chainx_original_magic_numbers, "In ChainX mode, use original magic numbers B = 100, alpha = 4 instead of variable B >= 100, alpha = 4")
		->group("")
		->check(CLI::Validator([&chain_p](std::string &_) {if (!chain_p.chainx and !chain_p.chainx_optimal) {return "pick --chainx or --chainx-opt";} else {return "";}}, "CHAINX(-OPT)", "ChainX mode"));

	app.add_flag("--chainx-opt-ensure-pred", chain_p.chainx_optimal_ensure_pred, "In ChainX-opt mode, always compute a ChainX-≺ chain at a tiny computational cost")
		->group("")
		->needs(chainx_optimal_option)
		->default_val(false);

	auto open_trunc_check = [](ofstream &out) {
		return CLI::Validator([&out](const string &file_path) {
			if (file_path.empty())
				return string("Empty path");
			out.open(file_path, std::ios::trunc);
			if (out.good())
				return string();
			else
				return string("Error in opening path: ") + file_path;
		}, "OUTFILE", "File can be opened and written into");
	};

	string mummer_output_path {""}, sam_output_path {""}, paf_output_path {""}, phylip_output_path {""};
	app.add_option("-o,--output", mummer_output_path, "Output the optimal chains in MUMmer-like format (ref start, query start, length) (1-based)")
		->check(open_trunc_check(ioss.mummer_out));
	auto sam_opt = app.add_option("-s,--sam", sam_output_path, "Output the optimal chains in SAM format")
		->check(open_trunc_check(ioss.sam_out));
	chain_p.store_sam_sequence = false;
	app.add_flag("--store-SAM-sequence", chain_p.store_sam_sequence, "Store the query sequence in the SAM output")
		->needs(sam_opt);
	app.add_option("-p,--paf", paf_output_path, "Output the optimal chains in PAF format")
		->check(open_trunc_check(ioss.paf_out));
	app.add_option("--phylip", phylip_output_path, "Output the anchored edit distances in PHYLIP format (lower triangular)")
		->check(open_trunc_check(ioss.phylip_out))
		->needs(all_to_all_option);

	long debug_random_anchors {0};
	auto debug_opt = app.add_option("--debug-random-anchors", debug_random_anchors, "DEBUG: fuzzy-test this number of random anchors in a 400x200 2D space")
		->group("")
		//->excludes(t_opt)
		//->excludes(q_opt)
		//->excludes(chainx_option)
		//->excludes(chainx_optimal_option)
		//->excludes(ata_opt)
		->check(CLI::Range(0L, std::numeric_limits<long>::max()));

	int debug_random_seed {-1};
	app.add_option("--debug-random-seed", debug_random_seed, "DEBUG: seed for random anchor generation (-1 is different at every invocation)")
		->group("")
		->default_val(-1)
		->needs(debug_opt);

	string debug_case_three_viz_path {""};
	app.add_option("--debug-case-three-viz", debug_case_three_viz_path, "DEBUG: visualize the space partitioning of case 3 in this output file (BMP)")
		->group("")
		->needs(debug_opt);

	try { app.parse(argc, argv); } catch(const CLI::ParseError &e) { return app.exit(e); }

	// ### additional checks ###
	if (debug_random_anchors > 0 and ((text_path != "") or (query_path != ""))) {
		return app.exit(CLI::Error("debug_mode", "--debug_random_anchors: excludes --text and --query", 1));
	}
	if (!all_to_all and debug_random_anchors == 0 and ((text_path == "") or (query_path == ""))) {
		return app.exit(CLI::Error("any_task", "Specify both a --text and a --query file", 1));
	}

	if (mode_s == "global") {
		chain_p.mode = utils::chaining_mode::global;
		chain_p.chainx_mode = chainx::mode::global;
	} else if (mode_s == "semiglobal") {
		chain_p.mode = utils::chaining_mode::semiglobal;
		chain_p.chainx_mode = chainx::mode::semiglobal;
	} else {
		return app.exit(CLI::Error("wrong_mode", "--mode: pick a correct chaining mode (global or semiglobal)", 1));
	}

	seed_p.anchor_type = MUM;
	if      (anchor_type_s == "MUM") seed_p.anchor_type = MUM;
	else if (anchor_type_s == "MEM") seed_p.anchor_type = MEM;
	else { return app.exit(CLI::Error("wrong_anchor_type", "--anchor-type: pick a supported anchor type (MUM or MEM)", 1)); };

	const bool backtrack = ioss.paf_out.is_open() or ioss.mummer_out.is_open() or ioss.sam_out.is_open();

	if ((text_path != "") or (query_path != "")) {
		// header for stdout output
		cout << "#text_id";
		cout << "\tquery_id";
		cout << "\tanchored_ED";
		cout << "\tanchors";
		cout << "\tmerged";
		cout << "\tseeding";
		cout << "\tpreprocessing";
		cout << "\tchaining";
		cout << (backtrack ? "\tbacktrack" : "");
		cout << "\ttotal_query";
		cout << ((chain_p.chainx or chain_p.chainx_optimal) ? "\trevisions" : "");
		cout << "\n";
	}

	// ### normal mode ###
	if ((text_path != "") and (query_path != "")) {
		cerr << string() +
			"[llchain] running with options" +
			" -t " + text_path +
			" -q " + query_path +
			((custom_anchors_path == "") ? string(" -m " + mode_s) : string("")) +
			((custom_anchors_path == "") ? string(" -a " + anchor_type_s) : string("")) +
			((custom_anchors_path == "") ? string(" -l " + to_string(seed_p.anchor_length)) : string("")) +
			((custom_anchors_path != "") ? string(" --custom-anchors " + custom_anchors_path) : string("")) +
			(reverse_complement ? string("-r") : string("")) +
			((mummer_output_path  != "") ? string(" -o " + mummer_output_path) : string("")) +
			((sam_output_path     != "") ? string(" -s " + sam_output_path) : string("")) +
			((paf_output_path     != "") ? string(" -p " + paf_output_path) : string(""))
			<< endl;

		htslib_wrapper::fasta_file text_idx(text_path);

		if (sam_output_path != "") {
			assert(ioss.sam_out.is_open());
			algo::write_SAM_header(ioss.sam_out);
			for (htslib_wrapper::hts_pos_t t = 0; t < text_idx.nseq(); t++)
				algo::write_SAM_text(text_idx.len(t), text_idx.id(t), ioss.sam_out);
		}

		if (mummer_output_path != "") {
			assert(ioss.mummer_out.is_open());
		}

		if (paf_output_path != "") {
			assert(ioss.paf_out.is_open());
		}

		ifstream custom_anchors_fs;
		if (custom_anchors_path != "") {
			custom_anchors_fs = ifstream(custom_anchors_path);
			auto found_custom_anchors_header = read_mummer_anchors_single(custom_anchors_fs).size();
			assert(found_custom_anchors_header == 0);
		}

		for (htslib_wrapper::hts_pos_t t = 0; t < text_idx.nseq(); t++) {
			auto start = std::chrono::steady_clock::now();
			const string text = text_idx.fetch(t);
			const string text_id = text_idx.id(t);
			auto const index = ((custom_anchors_path == "") ? mummer_essaMEM_wrapper::index(text, seed_p.anchor_length) : mummer_essaMEM_wrapper::dummy_index());
			const std::chrono::duration<double> index_time = std::chrono::steady_clock::now() - start;

			if (custom_anchors_path == "") {
				cerr << string() +
					"[llchain] indexed " + text_id +
					" in " + (std::ostringstream() << index_time).str()
					<< endl;
			}

			kseq::FastaGzInput fgz(query_path);
			string query_id, query;
			while (fgz.read_sequence(query_id, query)) {
				seed_and_chain(seed_p, chain_p, ioss, index, text, text_id, query, query_id, backtrack, false);

				if (reverse_complement) {
					llchain::utils::reverse_complement(query);
					seed_and_chain(seed_p, chain_p, ioss, index, text, text_id, query, query_id, backtrack, true);
				}
			}
		}

		if (ioss.sam_out.is_open()) {
			ioss.sam_out.close();
		}
		if (ioss.mummer_out.is_open()) {
			ioss.mummer_out.close();
		}
		if (ioss.paf_out.is_open()) {
			ioss.paf_out.close();
		}

		return 0;
	}

	// ### all-to-all mode ###
	if (all_to_all) {
		cerr << string() +
			"[llchain] running with options" +
			" --all-to-all" +
			" -t " + text_path +
			((custom_anchors_path == "") ? string(" -m " + mode_s) : string("")) +
			((custom_anchors_path == "") ? string(" -a " + anchor_type_s) : string("")) +
			((custom_anchors_path == "") ? string(" -l " + to_string(seed_p.anchor_length)) : string("")) +
			((custom_anchors_path != "") ? string(" --custom-anchors " + custom_anchors_path) : string("")) +
			(reverse_complement ? string("-r") : string("")) +
			((mummer_output_path  != "") ? string(" -o " + mummer_output_path) : string("")) +
			((sam_output_path     != "") ? string(" -s " + sam_output_path) : string("")) +
			((paf_output_path     != "") ? string(" -p " + paf_output_path) : string("")) +
			((phylip_output_path  != "") ? string(" --phylip " + phylip_output_path) : string(""))
			<< endl;

		htslib_wrapper::fasta_file text_idx(text_path);
		assert(text_idx.nseq() > 0);

		// https://phylipweb.github.io/phylip/doc/distance.html
		if (phylip_output_path != "") {
			ioss.phylip_out << text_idx.nseq() << "\n";
			if (text_idx.id(0).length() <= 10)
				ioss.phylip_out << text_idx.id(0) << string(10 - text_idx.id(0).length(), ' ');
			else
				ioss.phylip_out << text_idx.id(0).substr(0, 10);
			ioss.phylip_out << "\n";
		}

		for (htslib_wrapper::hts_pos_t t = 1; t < text_idx.nseq(); t++) {
			auto start = std::chrono::steady_clock::now();
			const string text = text_idx.fetch(t);
			const string text_id = text_idx.id(t);

			if (phylip_output_path != "") {
				if (text_id.length() <= 10)
					ioss.phylip_out << text_id << string(10 - text_id.length(), ' ');
				else
					ioss.phylip_out << text_id.substr(0, 10);
			}

			auto const index = ((custom_anchors_path == "") ? mummer_essaMEM_wrapper::index(text, seed_p.anchor_length) : mummer_essaMEM_wrapper::dummy_index());
			const std::chrono::duration<double> index_time = std::chrono::steady_clock::now() - start;

			if (custom_anchors_path == "") {
				cerr << string() +
					"[llchain] indexed " + text_id +
					" in " + (std::ostringstream() << index_time).str()
					<< endl;
			}

			ifstream custom_anchors_fs;
			if (custom_anchors_path != "") {
				custom_anchors_fs = ifstream(custom_anchors_path);
				auto found_custom_anchors_header = read_mummer_anchors_single(custom_anchors_fs).size();
				assert(found_custom_anchors_header == 0);
			}

			for (htslib_wrapper::hts_pos_t q = 0; q < t; q++) {
				string query = text_idx.fetch(q);
				const string query_id = text_idx.id(q);
				utils::anchor_index_t anchored_ed = seed_and_chain(seed_p, chain_p, ioss, index, text, text_id, query, query_id, backtrack, false);

				if (reverse_complement) {
					llchain::utils::reverse_complement(query);
					utils::anchor_index_t rc_anchored_ed = seed_and_chain(seed_p, chain_p, ioss, index, text, text_id, query, query_id, backtrack, false);
					anchored_ed = std::min(anchored_ed, rc_anchored_ed);
				}

				if (phylip_output_path != "") {
					ioss.phylip_out << " " << anchored_ed;
				}
			}
			if (phylip_output_path != "") {
				ioss.phylip_out << "\n";
			}
		}

		return 0;
	}

	// ### debug/fuzzy test mode ###
	if (debug_random_anchors > 0) {
		const int width = 400;
		const int height = 200;

		// anchor generation
		vector<anchor_t> anchors = random_anchors(width, height, debug_random_anchors, debug_random_seed);
		merge_perfect_chains(anchors); // only maximal anchors!
		place_dummy_anchors(width, height, anchors);
		weak_sort_anchors(anchors);

		// solve and plot case 3 recursions
		Image image(width, height);
		vector<long long> costs;
		vector<anchor_t> weak_chain;
		algo::weak_solve_naive(anchors, chain_p.mode, costs, weak_chain);
		plot_gap_gap_lower_diag(image, anchors, costs); // plot case 3 recursions
		plot_anchors(image, anchors); // plot all anchors
		plot_anchors(image, weak_chain, utils::defaults::selected_anchor_color); // recolor the optimal chain
		if (debug_case_three_viz_path != "") image.writeToFile(debug_case_three_viz_path);
		assert(costs.back() == algo::compute_chain_cost(weak_chain, chain_p.mode));

		// solve via the llchain algo and compare
		vector<long long> new_costs;
		algo::weak_solve_loglinear_debug(anchors, width, height, chain_p.mode, new_costs, costs);
		assert(new_costs.size() == costs.size() and new_costs.back() == costs.back());

		vector<anchor_t> chain;
		algo::weak_backtrack(anchors, costs, chain_p.mode, chain);
		cerr << "[llchain] DEBUG: optimal chain has cost     " << costs.back() << endl;
		cerr << "[llchain] DEBUG: backtracked chain has cost " << algo::compute_chain_cost(chain, chain_p.mode) << endl;
		assert(algo::compute_chain_cost(chain, chain_p.mode) == algo::compute_chain_cost(weak_chain, chain_p.mode));

		return 0;
	}

	return 1;
}

utils::anchor_index_t seed_and_chain(
	const seeding_params &seed_p,
	const chaining_params &chain_p,
	io_streams &ioss,
	const mummer_essaMEM_wrapper::sparseSA &index,
	const string &text,
	const string &text_id,
	string &query,
	const string &query_id,
	const bool backtrack,
	const bool is_query_rc
) {
	auto querystart = std::chrono::steady_clock::now();
	auto start = querystart;
	vector<anchor_t> matches;
	if (not ioss.custom_anchors_fs.is_open()) {
		if (seed_p.anchor_type == MUM)
			mummer_essaMEM_wrapper::find_MUMs(index, query, seed_p.anchor_length, matches);
		else if (seed_p.anchor_type == MEM)
			mummer_essaMEM_wrapper::find_MEMs(index, query, seed_p.anchor_length, matches);
		else
			assert(false);
	} else {
		matches = read_mummer_anchors_single(ioss.custom_anchors_fs);
	}
	const std::chrono::duration<double> seeding_time = std::chrono::steady_clock::now() - start;
	const long long found_anchors = matches.size();

	start = std::chrono::steady_clock::now();
	merge_perfect_chains(matches);
	place_dummy_anchors(text.length(), query.length(), matches);
	if (chain_p.chainx or chain_p.chainx_optimal) {
		chainx_sort_anchors(matches);
	} else {
		weak_sort_anchors(matches);
	}
	const std::chrono::duration<double> preprocessing_time = std::chrono::steady_clock::now() - start;

	start = std::chrono::steady_clock::now();
	vector<utils::anchor_index_t> costs;
	int chainx_revisions = -1;
	if (chain_p.chainx) {
		chainx::chainx(matches, query.length(), costs, chainx_revisions, chain_p.chainx_mode, false, chain_p.chainx_optimal_ensure_pred, chain_p.chainx_original_magic_numbers);
	} else if (chain_p.chainx_optimal) {
		chainx::chainx(matches, query.length(), costs, chainx_revisions, chain_p.chainx_mode, true,  chain_p.chainx_optimal_ensure_pred, chain_p.chainx_original_magic_numbers);
	} else {
		algo::weak_solve_loglinear(matches, text.length(), query.length(), chain_p.mode, costs);
	}
	const std::chrono::duration<double> main_chaining_time = std::chrono::steady_clock::now() - start;

	start = std::chrono::steady_clock::now();
	vector<anchor_t> chain;
	if (chain_p.chainx or (chain_p.chainx_optimal and chain_p.chainx_optimal_ensure_pred)) {
		algo::chainx_backtrack(matches, costs, chain_p.mode, chain);
	} else {
		algo::weak_backtrack(matches, costs, chain_p.mode, chain);
	}
	const std::chrono::duration<double> backtrack_chaining_time = std::chrono::steady_clock::now() - start;
	const std::chrono::duration<double> query_time = std::chrono::steady_clock::now() - querystart;

	// log chaining results
	cout << text_id << "\t" << query_id + (is_query_rc ? "_rc" : "");
	cout << "\t" << costs.back(); // anchored edit distance
	cout << "\t" << found_anchors;
	cout << "\t" << matches.size() - 2; // merged anchors
	cout << "\t" << seeding_time;
	cout << "\t" << preprocessing_time;
	cout << "\t" << main_chaining_time;
	if (backtrack) cout << "\t" << backtrack_chaining_time;
	cout << "\t" << query_time;
	cout << ((chain_p.chainx or chain_p.chainx_optimal) ? "\t" + to_string(chainx_revisions) : "");
	cout << "\n";

	if (ioss.sam_out.is_open() and chain.size() > 2) {
		algo::write_SAM_entry(text, text_id, query, query_id, is_query_rc, chain_p.store_sam_sequence, chain, chain_p.mode, costs.back(), ioss.sam_out);
	}
	if (ioss.paf_out.is_open() and chain.size() > 2) {
		algo::write_PAF_entry(text, text_id, query, query_id, is_query_rc, chain, chain_p.mode, costs.back(), ioss.paf_out);
	}
	if (ioss.mummer_out.is_open() and chain.size() > 2) {
		ioss.mummer_out << ">" << query_id + "_rc" << " (Reference " << text_id << ")\n";
		assert(chain.size() >= 2);
		for (size_type i = 1; i < chain.size() - 1; i++) {
			ioss.mummer_out << get<0>(chain[i])+1 << '\t' << get<1>(chain[i])+1 << '\t' << get<2>(chain[i])+1 << "\n";
		}
	}

	return costs.back();
}

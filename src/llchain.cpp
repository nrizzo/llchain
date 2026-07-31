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
import chainx;
import cli11;

using std::cout, std::cerr, std::endl;
using std::string, std::to_string;
using std::vector;
using std::ifstream, std::ofstream;
using namespace llchain;
using utils::anchor_t, utils::random_anchors, utils::plot_gap_gap_lower_diag, utils::plot_anchors, utils::Image, utils::place_dummy_anchors, utils::weak_sort_anchors, utils::chainx_sort_anchors, utils::merge_perfect_chains, utils::read_mummer_anchors_single;
typedef std::size_t size_type;

enum anchor_type_e { MUM, MEM };

int main(int argc, char **argv)
{
	// ### CLI setup ###
	CLI::App app("Log-linear chaining for anchored edit distance");
	string text_path {""}, query_path {""};
	auto t_opt = app.add_option("-t,--text", text_path, "One or more reference sequences in FASTA format (.fa,.fa.gz)")
		->check(CLI::ExistingFile);
	app.add_option("-q,--queries,--query", query_path, "Queries in FASTA format (.fa,.fa.gz)")
		->check(CLI::ExistingFile);

	string anchor_type_s {"MUM"};
	auto anchor_type_opt = app.add_option("-a,--anchor-type", anchor_type_s, "MUM or MEM (computed by MUMmer/essaMEM)");

	int anchor_length {20};
	auto anchor_length_opt = app.add_option("-l,--length", anchor_length, "Minimum anchor length")
		->check(CLI::Range(2, std::numeric_limits<int>::max()));

	string custom_anchors_path {""};
	app.add_option("--custom-anchors", custom_anchors_path, "Do not index/query but read the anchors from this file (NB it should respect the same order as query file)")
		->check(CLI::ExistingFile)
		->excludes(anchor_type_opt)
		->excludes(anchor_length_opt);

	string mode_s {""};
	app.add_option("-m,--mode", mode_s, "Chaining mode (global or semiglobal)")
		->default_val("global");

	bool all_to_all {false};
	auto ata_opt = app.add_flag("--all-to-all", all_to_all, "Pairwise comparisons (of the queries)")
		->excludes(t_opt);

	bool chainx {false};
	auto chainx_option = app.add_flag("--chainx", chainx, "Chain with at-cg/ChainX algorithm (variable B >= 100, alpha = 4)")
		->group("");

	bool chainx_optimal {false};
	auto chainx_optimal_option = app.add_flag("--chainx-opt", chainx_optimal, "Chain with algbio/ChainX *optimal* algorithm (variable B >= 100, alpha = 4)")
		->group("")
		->excludes(chainx_option);

	bool chainx_original_magic_numbers {false};
	app.add_flag("--chainx-original-magic-numbers", chainx_original_magic_numbers, "In ChainX mode, use original magic numbers B = 100, alpha = 4 instead of variable B >= 100, alpha = 4")
		->group("")
		->check(CLI::Validator([&chainx, &chainx_optimal](std::string &_) {if (!chainx and !chainx_optimal) {return "pick --chainx or --chainx-opt";} else {return "";}}, "CHAINX(-OPT)", "ChainX mode"));

	bool chainx_optimal_ensure_pred {false};
	app.add_flag("--chainx-opt-ensure-pred", chainx_optimal_ensure_pred, "In ChainX-opt mode, always compute a ChainX-≺ chain at a tiny computational cost")
		->group("")
		->needs(chainx_optimal_option);

	string mummer_output_path {""}, sam_output_path {""};
	app.add_option("-o,--output", mummer_output_path, "Output each optimal chain in MUMmer-like format {ref start, query start, length} (1-based)")
		->excludes(ata_opt);
	auto sam_opt = app.add_option("-s,--sam", sam_output_path, "Output approximate alignment based on the optimal chain (SAM format)")
		->excludes(ata_opt);

	bool store_sam_sequence {false};
	app.add_flag("--store-SAM-sequence", store_sam_sequence, "Store the query sequence in the SAM output")
		->needs(sam_opt);

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
		return app.exit(CLI::Error("any_task", "Specify both a text and a query file", 1));
	}

	algo::chaining_mode mode;
	chainx::mode chainx_mode;
	if (mode_s == "global") {
		mode = algo::chaining_mode::global;
		chainx_mode = chainx::mode::global;
	} else if (mode_s == "semiglobal") {
		mode = algo::chaining_mode::semiglobal;
		chainx_mode = chainx::mode::semiglobal;
	} else {
		return app.exit(CLI::Error("wrong_mode", "--mode: pick a correct chaining mode (global or semiglobal)", 1));
	}

	anchor_type_e anchortype;
	if      (anchor_type_s == "MUM") anchortype = MUM;
	else if (anchor_type_s == "MEM") anchortype = MEM;
	else { return app.exit(CLI::Error("wrong_anchor_type", "--anchor-type: pick a supported anchor type (MUM or MEM)", 1)); };

	cerr << "DEBUG: " << ((mode == algo::chaining_mode::global) ? "global" : "semiglobal") << " mode" << endl;
	if (custom_anchors_path == "") {
		cerr << "DEBUG: " << ((anchortype == MUM) ? "MUM" : "MEM") << " anchors of length >= " << anchor_length << endl;
	} else {
		cerr << "DEBUG: custom anchors from file " << custom_anchors_path << endl;
	}

	// ### normal mode ###
	if ((text_path != "") and (query_path != "")) {
		vector<string> texts, text_ids; // queries, query_ids;
		kseq::read_sequences(text_path, texts, text_ids);

		cerr << "DEBUG: read text sequences ";
		for (auto const &id : text_ids) cerr << id << " ";
		cerr << "of sizes ";
		for (auto const &text : texts) cerr << text.size() << " ";
		cerr << endl;

		ofstream sam_out;
		if (sam_output_path != "") {
			sam_out = ofstream(sam_output_path, std::ios::trunc);
			algo::write_SAM_header(sam_out);
			for (size_type t = 0; t < texts.size(); t++)
				algo::write_SAM_text(texts[t], text_ids[t], sam_out);
		}

		ofstream mummer_out;
		if (mummer_output_path != "") {
			mummer_out = ofstream(mummer_output_path, std::ios::trunc);
		}

		for (size_type t = 0; t < texts.size(); t++) {
			auto start = std::chrono::steady_clock::now(), querystart = std::chrono::steady_clock::now();
			auto index = ((custom_anchors_path == "") ? mummer_essaMEM_wrapper::index(texts[t], anchor_length) : mummer_essaMEM_wrapper::dummy_index());
			const std::chrono::duration<double> index_time = std::chrono::steady_clock::now() - start;
			if (custom_anchors_path == "") {
				cerr << "DEBUG: indexed " << text_ids[t] << " in " << index_time << endl;
			}
			ifstream custom_anchors_fs;
			if (custom_anchors_path != "") {
				custom_anchors_fs = ifstream(custom_anchors_path);
				auto a = read_mummer_anchors_single(custom_anchors_fs);
				assert(a.size() == 0);
			}

			kseq::FastaGzInput fgz(query_path);
			string query_id, query;
			while (fgz.read_sequence(query_id, query)) {
				cerr << "DEBUG: querying " << query_id << " in " << text_ids[t] << " (" << ((anchortype == MUM) ? "MUM" : "MEM") << " seeds of length >= " << anchor_length << ")...";

				querystart = std::chrono::steady_clock::now();
				start = querystart;
				vector<anchor_t> matches;
				if (custom_anchors_path  == "") {
					if (anchortype == MUM)
						mummer_essaMEM_wrapper::find_MUMs(index, query, anchor_length, matches);
					else if (anchortype == MEM)
						mummer_essaMEM_wrapper::find_MEMs(index, query, anchor_length, matches);
				} else {
					matches = read_mummer_anchors_single(custom_anchors_fs);
				}
				const std::chrono::duration<double> seeding_time = std::chrono::steady_clock::now() - start;
				const long long found_anchors = matches.size();

				start = std::chrono::steady_clock::now();
				merge_perfect_chains(matches);
				place_dummy_anchors(texts[t].size(), query.size(), matches);
				if (chainx or chainx_optimal) {
					chainx_sort_anchors(matches);
				} else {
					weak_sort_anchors(matches);
				}
				const std::chrono::duration<double> preprocessing_time = std::chrono::steady_clock::now() - start;

				start = std::chrono::steady_clock::now();
				vector<utils::anchor_index_t> costs;
				int chainx_revisions = -1;
				if (chainx) {
					chainx::chainx(matches, query.size(), costs, chainx_revisions, chainx_mode, false, chainx_optimal_ensure_pred, chainx_original_magic_numbers);
				} else if (chainx_optimal) {
					chainx::chainx(matches, query.size(), costs, chainx_revisions, chainx_mode, true,  chainx_optimal_ensure_pred, chainx_original_magic_numbers);
				} else {
					algo::weak_solve_loglinear(matches, texts[t].size(), query.size(), mode, costs);
				}
				const std::chrono::duration<double> main_chaining_time = std::chrono::steady_clock::now() - start;

				start = std::chrono::steady_clock::now();
				vector<anchor_t> chain;
				if (chainx or (chainx_optimal and chainx_optimal_ensure_pred)) {
					algo::chainx_backtrack(matches, costs, mode, chain);
				} else {
					algo::weak_backtrack(matches, costs, mode, chain);
				}
				const std::chrono::duration<double> backtrack_chaining_time = std::chrono::steady_clock::now() - start;
				const std::chrono::duration<double> query_time = std::chrono::steady_clock::now() - querystart;

				cerr << "done (" << found_anchors << " anchors, " << matches.size() - 2 << " merged, " << costs.back() << " anchored edit distance, " << seeding_time << " seeding, " << preprocessing_time << " preprocessing, " << main_chaining_time << " chaining, " << backtrack_chaining_time << " backtrack, " << query_time << " total query time" << ((chainx or chainx_optimal) ? (", " + to_string(chainx_revisions) + " revisions") : "") << ")" << endl;

				if (sam_output_path != "" and chain.size() > 2) {
					algo::write_SAM_entry(texts[t], text_ids[t], query, query_id, store_sam_sequence, chain, mode, costs.back(), sam_out);
				}
				if (mummer_output_path != "") {
					mummer_out << ">" << query_id << " (Reference " << text_ids[t] << ")\n";
					assert(chain.size() >= 2);
					for (size_type i = 1; i < chain.size() - 1; i++) {
						mummer_out << get<0>(chain[i])+1 << '\t' << get<1>(chain[i])+1 << '\t' << get<2>(chain[i])+1 << "\n";
					}
				}
			}
		}
		if (sam_output_path != "") {
			sam_out.close();
		}
		if (mummer_output_path != "") {
			mummer_out.close();
		}

		return 0;
	}

	// ### all-to-all mode ###
	if (all_to_all) {
		vector<string> queries, query_ids;
		kseq::read_sequences(query_path, queries, query_ids);

		vector<vector<utils::anchor_index_t>> distances(
				queries.size(),
				vector<utils::anchor_index_t>(queries.size(),
					std::numeric_limits<utils::anchor_index_t>::max())
				);
		for (size_type i = 1; i < queries.size(); i++) {
			auto start = std::chrono::steady_clock::now(), querystart = std::chrono::steady_clock::now();
			auto const index = mummer_essaMEM_wrapper::index(queries[i], anchor_length);
			const std::chrono::duration<double> index_time = std::chrono::steady_clock::now() - start;
			cerr << "DEBUG: indexed " << query_ids[i] << " in " << index_time << endl;

			for (size_type j = 0; j < i; j++) {
				cerr << "DEBUG: querying " << query_ids[j] << " in " << query_ids[i] << " (" << ((anchortype == MUM) ? "MUM" : "MEM") << " seeds of length >= " << anchor_length << ")...";
				querystart = std::chrono::steady_clock::now();
				start = querystart;
				vector<anchor_t> matches;
				if (anchortype == MUM)
					mummer_essaMEM_wrapper::find_MUMs(index, queries[j], anchor_length, matches);
				else if (anchortype == MEM)
					mummer_essaMEM_wrapper::find_MEMs(index, queries[j], anchor_length, matches);
				const std::chrono::duration<double> seeding_time = std::chrono::steady_clock::now() - start;
				const long long found_anchors = matches.size();

				start = std::chrono::steady_clock::now();
				merge_perfect_chains(matches);
				place_dummy_anchors(queries[i].size(), queries[j].size(), matches);
				if (chainx or chainx_optimal) {
					chainx_sort_anchors(matches);
				} else {
					weak_sort_anchors(matches);
				}
				const std::chrono::duration<double> preprocessing_time = std::chrono::steady_clock::now() - start;

				start = std::chrono::steady_clock::now();
				vector<utils::anchor_index_t> costs;
				int chainx_revisions = -1;
				if (chainx) {
					chainx::chainx(matches, queries[j].size(), costs, chainx_revisions, chainx_mode, false);
				} else if (chainx_optimal) {
					chainx::chainx(matches, queries[j].size(), costs, chainx_revisions, chainx_mode);
				} else {
					algo::weak_solve_loglinear(matches, queries[i].size(), queries[j].size(), mode, costs);
				}
				const std::chrono::duration<double> main_chaining_time = std::chrono::steady_clock::now() - start;
				const std::chrono::duration<double> query_time = std::chrono::steady_clock::now() - querystart;
				distances[i][j] = costs.back();

				cerr << "done (" << found_anchors << " anchors, " << matches.size() - 2 << " merged, " << costs.back() << " anchored edit distance, " << seeding_time << " seeding, " << preprocessing_time << " preprocessing, " << main_chaining_time << " chaining, " << query_time << " total query time" << ((chainx or chainx_optimal) ? (", " + to_string(chainx_revisions) + " revisions") : "") << ")" << endl;
			}
		}

		// https://phylipweb.github.io/phylip/doc/distance.html
		cout << queries.size() << "\n";
		for (size_type i = 0; i < queries.size(); i++) {
			if (query_ids[i].size() <= 10)
				cout << query_ids[i] << string(10 - query_ids[i].size(), ' ');
			else
				cout << query_ids[i].substr(0, 10);

			for (size_type j = 0; j < queries.size(); j++) {
				if (i == j)
					cout << " " << 0;
				else
					cout << " " << ((j < i) ? distances[i][j] : distances[j][i]);
			}

			cout << "\n";
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
		algo::weak_solve_naive(anchors, mode, costs, weak_chain);
		plot_gap_gap_lower_diag(image, anchors, costs); // plot case 3 recursions
		plot_anchors(image, anchors); // plot all anchors
		plot_anchors(image, weak_chain, utils::defaults::selected_anchor_color); // recolor the optimal chain
		image.writeToFile(debug_case_three_viz_path);
		assert(costs.back() == algo::compute_chain_cost(weak_chain, mode));

		// solve via the llchain algo and compare
		vector<long long> new_costs;
		algo::weak_solve_loglinear_debug(anchors, width, height, mode, new_costs, costs);
		assert(new_costs.size() == costs.size() and new_costs.back() == costs.back());

		vector<anchor_t> chain;
		algo::weak_backtrack(anchors, costs, mode, chain);
		cerr << "Optimal chain has cost     " << costs.back() << endl;
		cerr << "Backtracked chain has cost " << algo::compute_chain_cost(chain, mode) << endl;
		assert(algo::compute_chain_cost(chain, mode) == algo::compute_chain_cost(weak_chain, mode));

		return 0;
	}

	return 1;
}

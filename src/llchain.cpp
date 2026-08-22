#ifndef VERSION
#define VERSION "GitHub"
#endif

module;
#include <chrono> // fix for GCC 15.2.1
#include "concurrentqueue.h"

import <iostream>;
import <string>;
import <sstream>;
import <vector>;
import <fstream>;
#pragma GCC diagnostic push // GCC 16.1.1
#pragma GCC diagnostic ignored "-Wsign-compare"
import <syncstream>;
#pragma GCC diagnostic pop
import <thread>;
import <latch>;
import <barrier>;
import <cassert>;
import <optional>;
import <atomic>;
import <algorithm>;
import utils;
import algo;
import kseq;
import mummer_essaMEM_wrapper;
import htslib_wrapper;
import chainx;
import cli11;

using std::cout, std::cerr, std::endl;
using std::string, std::to_string;
using std::ostringstream;
using std::vector;
using std::ifstream, std::ofstream;
using std::osyncstream;
using std::jthread;
using std::latch;
using std::barrier;
using std::optional;
using std::atomic;
using std::min;
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
	ofstream sam_out;
	ofstream paf_out;
	ofstream mummer_out;
};
struct io_streams_sync {
	osyncstream cout;
	osyncstream sam_out;
	osyncstream paf_out;
	osyncstream mummer_out;
};
struct read {
	string read;
	string read_id;
};
struct ata_read {
	string read;
	string read_id;
	htslib_wrapper::hts_pos_t read_i;
};
struct assembly_traits : public moodycamel::ConcurrentQueueDefaultTraits {
	static const size_t BLOCK_SIZE = 2;
};

/*
 * seed and compute the anchored edit distance
 */
utils::anchor_index_t seed_and_chain(
	const seeding_params &seed_p,
	const chaining_params &chain_p,
	io_streams_sync &ioss,
	const mummer_essaMEM_wrapper::sparseSA &index,
	const string::size_type text_length,
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
	int qthreads;
	app.add_option("--qthreads", qthreads, "Query threads")
		->check(CLI::PositiveNumber)
		->default_val(1);


	string anchor_type_s {}; // to be replaced by seed_p.anchor_type
	app.add_option("-a,--anchor-type", anchor_type_s, "MUM or MEM (computed via MUMmer/essaMEM)")
		->default_val("MUM");

	app.add_option("-l,--length", seed_p.anchor_length, "Minimum anchor length")
		->check(CLI::Range(2, std::numeric_limits<int>::max()))
		->default_val(20);

	//string custom_anchors_path {""};
	//app.add_option("--custom-anchors", custom_anchors_path, "Do not index/query but read the anchors from this file (NB it should respect the same order as query file)")
	//	->check(CLI::ExistingFile)
	//	->excludes(anchor_type_opt)
	//	->excludes(anchor_length_opt);

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
		->check(CLI::Validator([&chain_p](string &_) {if (!chain_p.chainx and !chain_p.chainx_optimal) {return "pick --chainx or --chainx-opt";} else {return "";}}, "CHAINX(-OPT)", "ChainX mode"));

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
	ofstream phylip_out;
	app.add_option("--phylip", phylip_output_path, "Output the distances in a PHYLIP matrix (lower triangular)")
		->check(open_trunc_check(phylip_out))
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
			" --qthreads " + to_string(qthreads) +
			string(" -m " + mode_s) +
			string(" -a " + anchor_type_s) +
			string(" -l " + to_string(seed_p.anchor_length)) +
			(reverse_complement ? string("-r") : string("")) +
			((mummer_output_path  != "") ? string(" -o " + mummer_output_path) : string("")) +
			((sam_output_path     != "") ? string(" -s " + sam_output_path) : string("")) +
			((paf_output_path     != "") ? string(" -p " + paf_output_path) : string(""))
			<< endl;

		htslib_wrapper::fasta_file text_idx(text_path);

		if (mummer_output_path != "") {
			assert(ioss.mummer_out.is_open());
		}

		if (paf_output_path != "") {
			assert(ioss.paf_out.is_open());
		}

		if (sam_output_path != "") {
			assert(ioss.sam_out.is_open());
			algo::write_SAM_header(ioss.sam_out);
			for (htslib_wrapper::hts_pos_t t = 0; t < text_idx.nseq(); t++)
				algo::write_SAM_text(text_idx.len(t), text_idx.id(t), ioss.sam_out);
		}

		// index all text sequences
		vector<string> texts(text_idx.nseq());
		vector<mummer_essaMEM_wrapper::sparseSA> text_indexes;
		text_indexes.reserve(text_idx.nseq());
		for (htslib_wrapper::hts_pos_t t = 0; t < text_idx.nseq(); t++) {
			auto start = std::chrono::steady_clock::now();
			texts[t] = text_idx.fetch(t);
			const string text_id = text_idx.id(t);
			text_indexes.push_back(mummer_essaMEM_wrapper::index(texts[t], seed_p.anchor_length));
			const std::chrono::duration<double> index_time = std::chrono::steady_clock::now() - start;

			cerr << string() +
				"[llchain] indexed " + text_id +
				" in " + (ostringstream() << index_time).str()
				<< endl;
		}

		latch reader_done(1);
		// queue for normal/short sequences
		moodycamel::ConcurrentQueue<read> read_queue(3*qthreads, 1, 1);
		// queue for very long sequences
		moodycamel::ConcurrentQueue<read,assembly_traits> assembly_queue(1, 1, 1);
		auto reader_thread = jthread([&query_path, &read_queue, &assembly_queue, &reader_done]() {
			kseq::FastaGzInput fgz(query_path);
			string query_id, query;
			while (fgz.read_sequence(query_id, query)) {
				assert(query.length() > 0 and query_id != "");
				if (query.length() <= 200000) {
					while(!read_queue.try_enqueue({query, query_id}))
						std::this_thread::sleep_for(std::chrono::milliseconds(20));
				} else {
					while(!assembly_queue.try_enqueue({query, query_id}))
						std::this_thread::sleep_for(std::chrono::milliseconds(20));
				}
			}
			reader_done.count_down();
		});

		vector<jthread> seed_chain_workers(qthreads);
		for (int i = 0; i < qthreads; i++) {
			seed_chain_workers[i] = jthread([&]() {
				read r;
				bool final_check = false;
				io_streams_sync ioss_sync{
					osyncstream(cout),
					osyncstream(ioss.sam_out),
					osyncstream(ioss.paf_out),
					osyncstream(ioss.mummer_out)
				};
				while (true) {
					while (read_queue.try_dequeue(r) or assembly_queue.try_dequeue(r)) {
						assert(r.read.length() > 0 and r.read_id != "");
						for (vector<string::size_type>::size_type t = 0; t < texts.size(); t++) {
							seed_and_chain(seed_p, chain_p, ioss_sync, text_indexes[t], texts[t].length(), text_idx.id(t), r.read, r.read_id, backtrack, false);

							if (reverse_complement) {
								llchain::utils::reverse_complement(r.read);
								seed_and_chain(seed_p, chain_p, ioss_sync, text_indexes[t], texts[t].length(), text_idx.id(t), r.read, r.read_id, backtrack, true);
							}
						}
						if (ioss_sync.cout.good()) ioss_sync.cout.emit();
						if (ioss_sync.sam_out.good()) ioss_sync.sam_out.emit();
						if (ioss_sync.paf_out.good()) ioss_sync.paf_out.emit();
						if (ioss_sync.mummer_out.good()) ioss_sync.mummer_out.emit();
					}
					if (reader_done.try_wait()) {
						if (final_check)
							break;
						final_check = true;
					}
					if (!final_check)
						std::this_thread::sleep_for(std::chrono::milliseconds(20));
				}
			});
		}

		reader_thread.join();
		for (auto &t : seed_chain_workers) {
			t.join();
		}

		return 0;
	}

	// ### all-to-all mode ###
	if (all_to_all) {
		cerr << string() +
			"[llchain] running with options" +
			" --all-to-all" +
			" -t " + text_path +
			" --qthreads " + to_string(qthreads) +
			string(" -m " + mode_s) +
			string(" -a " + anchor_type_s) +
			string(" -l " + to_string(seed_p.anchor_length)) +
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
			phylip_out << text_idx.nseq() << "\n";
			if (text_idx.id(0).length() <= 10)
				phylip_out << text_idx.id(0) << string(10 - text_idx.id(0).length(), ' ');
			else
				phylip_out << text_idx.id(0).substr(0, 10);
			phylip_out << "\n";
		}

		// queue for normal/short sequences
		moodycamel::ConcurrentQueue<ata_read> ata_read_queue(3*qthreads, 1, 1);
		// queue for very long sequences
		moodycamel::ConcurrentQueue<ata_read,assembly_traits> ata_assembly_queue(1, 1, 1);
		vector<utils::anchor_index_t> distances;
		barrier sync_point(1 + qthreads); // to sync reader_indexer and seed_chain workers
		atomic<bool> stage_end(false);
		string text(""), text_id("");
		optional<mummer_essaMEM_wrapper::sparseSA> index = mummer_essaMEM_wrapper::dummy_index(); // no copy/move constructor
		string::size_type text_length = 0;
		latch reader_done(1);
		auto reader_indexer_worker = jthread([&]() {
			for (htslib_wrapper::hts_pos_t t = 1; t < text_idx.nseq(); t++) {
				auto start = std::chrono::steady_clock::now();
				text = text_idx.fetch(t);
				text_id = text_idx.id(t);

				if (phylip_output_path != "") {
					if (text_id.length() <= 10)
						phylip_out << text_id << string(10 - text_id.length(), ' ');
					else
					phylip_out << text_id.substr(0, 10);
				}

				index.reset();
				index.emplace(mummer_essaMEM_wrapper::index(text, seed_p.anchor_length));
				text_length = text.length();
				const std::chrono::duration<double> index_time = std::chrono::steady_clock::now() - start;

				cerr << string() +
					"[llchain] indexed " + text_id +
					" in " + (ostringstream() << index_time).str()
					<< endl;

				if (phylip_output_path != "") {
					distances.resize(t);
				}
				for (htslib_wrapper::hts_pos_t q = 0; q < t; q++) {
					string query = text_idx.fetch(q);
					const string query_id = text_idx.id(q);
					if (query.length() <= 200000) {
						while (!ata_read_queue.try_enqueue({query, query_id, q}))
							std::this_thread::sleep_for(std::chrono::milliseconds(20));
					} else {
						while (!ata_assembly_queue.try_enqueue({query, query_id, q}))
							std::this_thread::sleep_for(std::chrono::milliseconds(20));
					}
				}

				stage_end = true;
				sync_point.arrive_and_wait(); // wait for the seed_chain_workers
				stage_end = false;

				if (phylip_output_path != "") {
					for (auto d : distances)
						phylip_out << " " << d;
					phylip_out << "\n";
				}
			}
			reader_done.count_down();
		});

		vector<jthread> seed_chain_workers(qthreads);
		for (int i = 0; i < qthreads; i++) {
			seed_chain_workers[i] = jthread([&]() {
				ata_read r;
				bool final_check = false;
				io_streams_sync ioss_sync{
					osyncstream(cout),
					osyncstream(ioss.sam_out),
					osyncstream(ioss.paf_out),
					osyncstream(ioss.mummer_out)
				};
				while (true) {
					while (ata_read_queue.try_dequeue(r) or ata_assembly_queue.try_dequeue(r)) {
						if (chain_p.mode == utils::chaining_mode::global or text_length >= r.read.length()) {
							auto anchored_ed = seed_and_chain(seed_p, chain_p, ioss_sync, *index, text_length, text_id, r.read, r.read_id, backtrack, false);
							if (reverse_complement) {
								llchain::utils::reverse_complement(r.read);
								anchored_ed = min(anchored_ed, seed_and_chain(seed_p, chain_p, ioss_sync, *index, text_length, text_id, r.read, r.read_id, backtrack, true));
							}
							if (phylip_output_path != "") {
								distances[r.read_i] = anchored_ed;
							}
						} else {
							// the smaller text sequence is the "query"
							auto qindex = mummer_essaMEM_wrapper::index(r.read, seed_p.anchor_length);
							auto anchored_ed = seed_and_chain(seed_p, chain_p, ioss_sync, qindex, r.read.length(), r.read_id, text, text_id, backtrack, false);
							if (reverse_complement) {
								string rc_text(text);
								llchain::utils::reverse_complement(rc_text);
								anchored_ed = min(anchored_ed, seed_and_chain(seed_p, chain_p, ioss_sync, qindex, r.read.length(), r.read_id, rc_text, text_id, backtrack, true));
							}
							if (phylip_output_path != "") {
								distances[r.read_i] = anchored_ed;
							}
						}
						if (ioss_sync.cout.good()) ioss_sync.cout.emit();
						if (ioss_sync.sam_out.good()) ioss_sync.sam_out.emit();
						if (ioss_sync.paf_out.good()) ioss_sync.paf_out.emit();
						if (ioss_sync.mummer_out.good()) ioss_sync.mummer_out.emit();
					}
					if (final_check and stage_end) {
						sync_point.arrive_and_wait();
						final_check = false;
						stage_end = false;
						continue;
					}
					if (reader_done.try_wait() or stage_end) {
						if (final_check and reader_done.try_wait())
							break;
						final_check = true;
					}
					if (!final_check)
						std::this_thread::sleep_for(std::chrono::milliseconds(20));
				}
			});
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
	io_streams_sync &ioss,
	const mummer_essaMEM_wrapper::sparseSA &index,
	const string::size_type text_length,
	const string &text_id,
	string &query,
	const string &query_id,
	const bool backtrack,
	const bool is_query_rc
) {
	auto querystart = std::chrono::steady_clock::now();
	auto start = querystart;
	vector<anchor_t> matches;
	if (seed_p.anchor_type == MUM)
		mummer_essaMEM_wrapper::find_MUMs(index, query, seed_p.anchor_length, matches);
	else if (seed_p.anchor_type == MEM)
		mummer_essaMEM_wrapper::find_MEMs(index, query, seed_p.anchor_length, matches);
	else
		assert(false);
	const std::chrono::duration<double> seeding_time = std::chrono::steady_clock::now() - start;
	const long long found_anchors = matches.size();

	start = std::chrono::steady_clock::now();
	merge_perfect_chains(matches);
	place_dummy_anchors(text_length, query.length(), matches);
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
		algo::weak_solve_loglinear(matches, text_length, query.length(), chain_p.mode, costs);
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
	ioss.cout << text_id << "\t" << query_id + (is_query_rc ? "_rc" : "");
	ioss.cout << "\t" << costs.back(); // anchored edit distance
	ioss.cout << "\t" << found_anchors;
	ioss.cout << "\t" << matches.size() - 2; // merged anchors
	ioss.cout << "\t" << seeding_time;
	ioss.cout << "\t" << preprocessing_time;
	ioss.cout << "\t" << main_chaining_time;
	if (backtrack) ioss.cout << "\t" << backtrack_chaining_time;
	ioss.cout << "\t" << query_time;
	ioss.cout << ((chain_p.chainx or chain_p.chainx_optimal) ? "\t" + to_string(chainx_revisions) : "");
	ioss.cout << "\n";

	if (ioss.sam_out.good() and chain.size() > 2) {
		algo::write_SAM_entry(text_length, text_id, query, query_id, is_query_rc, chain_p.store_sam_sequence, chain, chain_p.mode, costs.back(), ioss.sam_out);
	}
	if (ioss.paf_out.good() and chain.size() > 2) {
		algo::write_PAF_entry(text_length, text_id, query, query_id, is_query_rc, chain, chain_p.mode, costs.back(), ioss.paf_out);
	}
	if (ioss.mummer_out.good() and chain.size() > 2) {
		ioss.mummer_out << ">" << query_id + "_rc" << " (Reference " << text_id << ")\n";
		assert(chain.size() >= 2);
		for (size_type i = 1; i < chain.size() - 1; i++) {
			ioss.mummer_out << get<0>(chain[i])+1 << '\t' << get<1>(chain[i])+1 << '\t' << get<2>(chain[i])+1 << "\n";
		}
	}

	return costs.back();
}

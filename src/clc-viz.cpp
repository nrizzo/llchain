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

#include "command-line-parsing/cmdline.h" // cmdline_parser (gengetopt)

using std::cout, std::cerr, std::endl;
using std::string;
using std::vector;
using std::ifstream;
using utils::anchor_t, utils::random_anchors, utils::plot_gap_gap_lower_diag, utils::plot_anchors, utils::Image, utils::place_dummy_anchors, utils::sort_anchors, utils::merge_perfect_chains, utils::read_mummer_anchors_single;
typedef std::size_t size_type;

enum anchor_type { MUM, MEM };

int main(int argc, char **argv) 
{
	// setup
	gengetopt_args_info argsinfo;
	if (cmdline_parser(argc, argv, &argsinfo) != 0) exit(1);
	if (argsinfo.text_arg == NULL and argsinfo.query_arg == NULL and argsinfo.random_anchors_arg == -1) {
		cmdline_parser_print_help();
		exit(1);
	}

	// check parameters
	if (!argsinfo.all_to_all_flag and (argsinfo.text_arg == NULL) xor (argsinfo.query_arg == NULL)) {
		cerr << "Error: in normal mode, specify both a text and a query file!" << endl;
		exit(1);
	}
	if (((argsinfo.text_arg != NULL) or (argsinfo.query_arg != NULL)) and argsinfo.random_anchors_arg != -1) {
		cerr << "Error: input files and random anchor generation are not compatible!" << endl;
		exit(1);
	}
	if (argsinfo.all_to_all_flag and ((argsinfo.text_arg != NULL) or (argsinfo.query_arg == NULL) or (argsinfo.random_anchors_arg > 0))) {
		cerr << "Error: pairwise comparison mode (--all-to-all) accepts option --query only!" << endl;
		exit(1);
	}
	if (argsinfo.all_to_all_flag and argsinfo.custom_anchors_arg != NULL) {
		cerr << "Error: --custom-anchors is not compatible with --all-to-all!" << endl;
		exit(1);
	}

	const int anchorlength = argsinfo.anchor_length_arg;
	if (anchorlength <= 0) { cerr << "Error: pick an anchor length >= 1." << endl; exit(1); };

	if ((argsinfo.text_arg != NULL) and (argsinfo.query_arg != NULL)) {
		algo::chaining_mode mode;
		if      (string(argsinfo.mode_arg) == "global")     mode = algo::chaining_mode::global;
		else if (string(argsinfo.mode_arg) == "semiglobal") mode = algo::chaining_mode::semiglobal;
		else { cerr << "Error: pick a correct chaining mode (global/semiglobal)." << endl; exit(1); }

		anchor_type anchortype;
		if      (string(argsinfo.anchor_type_arg) == "MUM") anchortype = MUM;
		else if (string(argsinfo.anchor_type_arg) == "MEM") anchortype = MEM;
		else { cerr << "Error: pick a correct anchor type (MUM/MEM)." << endl; exit(1); }

		vector<string> texts, text_ids; // queries, query_ids;
		kseq::read_sequences(string(argsinfo.text_arg),  texts,   text_ids);

		cerr << "DEBUG: read text sequences ";
		for (auto const &id : text_ids) cerr << id << " ";
		cerr << "of sizes ";
		for (auto const &text : texts) cerr << text.size() << " ";
		cerr << endl;

		for (size_type t = 0; t < texts.size(); t++) {
			auto start = std::chrono::steady_clock::now(), querystart = std::chrono::steady_clock::now();
			auto index = ((argsinfo.custom_anchors_arg == NULL) ? mummer_essaMEM_wrapper::index(texts[t], anchorlength) : mummer_essaMEM_wrapper::dummy_index());
			const std::chrono::duration<double> index_time = std::chrono::steady_clock::now() - start;
			if (argsinfo.custom_anchors_arg == NULL) {
				cerr << "DEBUG: indexed " << text_ids[t] << " in " << index_time << endl;
			}
			ifstream custom_anchors_fs;
			if (argsinfo.custom_anchors_arg != NULL) {
				custom_anchors_fs = ifstream(string(argsinfo.custom_anchors_arg));
				auto a = read_mummer_anchors_single(custom_anchors_fs);
				assert(a.size() == 0);
			}

			kseq::FastaGzInput fgz(string(argsinfo.query_arg));
			string query_id, query;
			while (fgz.read_sequence(query_id, query)) {
				cerr << "DEBUG: querying " << query_id << " in " << text_ids[t] << " (" << ((anchortype == MUM) ? "MUM" : "MEM") << " seeds of length >= " << anchorlength << ")...";

				querystart = std::chrono::steady_clock::now();
				start = querystart;
				vector<anchor_t> matches;
				if (argsinfo.custom_anchors_arg == NULL) {
					if (anchortype == MUM)
						mummer_essaMEM_wrapper::find_MUMs(index, query, anchorlength, matches);
					else if (anchortype == MEM)
						mummer_essaMEM_wrapper::find_MEMs(index, query, anchorlength, matches);
				} else {
					matches = read_mummer_anchors_single(custom_anchors_fs);
				}
				const std::chrono::duration<double> seeding_time = std::chrono::steady_clock::now() - start;
				const long long found_anchors = matches.size();

				start = std::chrono::steady_clock::now();
				merge_perfect_chains(matches);
				place_dummy_anchors(texts[t].size(), query.size(), matches);
				sort_anchors(matches);
				const std::chrono::duration<double> preprocessing_time = std::chrono::steady_clock::now() - start;

				start = std::chrono::steady_clock::now();
				vector<utils::anchor_index_t> costs;
				algo::solve_linearithmic(matches, texts[t].size(), query.size(), mode, costs);
				vector<anchor_t> chain;
				algo::weak_backtrack(matches, costs, mode, chain);
				const std::chrono::duration<double> chaining_time = std::chrono::steady_clock::now() - start;
				const std::chrono::duration<double> query_time = std::chrono::steady_clock::now() - querystart;

				cerr << "done (" << found_anchors << " anchors, " << matches.size() - 2 << " merged, " << costs.back() << " anchored edit distance, " << seeding_time << " seeding, " << preprocessing_time << " preprocessing, " << chaining_time << " chaining, " << query_time << " total query time)" << endl;
			}
		}

		return 0;
	}

	if (argsinfo.all_to_all_flag) {
		algo::chaining_mode mode;
		if      (string(argsinfo.mode_arg) == "global")     mode = algo::chaining_mode::global;
		else if (string(argsinfo.mode_arg) == "semiglobal") mode = algo::chaining_mode::semiglobal;
		else { cerr << "Error: pick a correct chaining mode (global/semiglobal)." << endl; exit(1); }

		anchor_type anchortype;
		if      (string(argsinfo.anchor_type_arg) == "MUM") anchortype = MUM;
		else if (string(argsinfo.anchor_type_arg) == "MEM") anchortype = MEM;
		else { cerr << "Error: pick a correct anchor type (MUM/MEM)." << endl; exit(1); }

		vector<string> queries, query_ids;
		kseq::read_sequences(string(argsinfo.query_arg), queries, query_ids);

		vector<vector<utils::anchor_index_t>> distances(
				queries.size(),
				vector<utils::anchor_index_t>(queries.size(),
					std::numeric_limits<utils::anchor_index_t>::max())
				);
		for (size_type i = 0; i < queries.size(); i++) {
			auto start = std::chrono::steady_clock::now(), querystart = std::chrono::steady_clock::now();
			auto const index = mummer_essaMEM_wrapper::index(queries[i], anchorlength);
			const std::chrono::duration<double> index_time = std::chrono::steady_clock::now() - start;
			cerr << "DEBUG: indexed " << query_ids[i] << " in " << index_time << endl;

			for (size_type j = 0; j < i; j++) {
				querystart = std::chrono::steady_clock::now();
				start = querystart;
				vector<anchor_t> matches;
				if (anchortype == MUM)
					mummer_essaMEM_wrapper::find_MUMs(index, queries[j], anchorlength, matches);
				else if (anchortype == MEM)
					mummer_essaMEM_wrapper::find_MEMs(index, queries[j], anchorlength, matches);
				const std::chrono::duration<double> seeding_time = std::chrono::steady_clock::now() - start;
				const long long found_anchors = matches.size();

				start = std::chrono::steady_clock::now();
				merge_perfect_chains(matches);
				place_dummy_anchors(queries[i].size(), queries[j].size(), matches);
				sort_anchors(matches);
				const std::chrono::duration<double> preprocessing_time = std::chrono::steady_clock::now() - start;

				start = std::chrono::steady_clock::now();
				vector<utils::anchor_index_t> costs;
				algo::solve_linearithmic(matches, queries[i].size(), queries[j].size(), mode, costs);
				const std::chrono::duration<double> chaining_time = std::chrono::steady_clock::now() - start;
				const std::chrono::duration<double> query_time = std::chrono::steady_clock::now() - querystart;
				distances[i][j] = costs.back();
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

	if (argsinfo.random_anchors_arg > 0) {
		const int width = 400;
		const int height = 200;

		// anchor generation
		vector<anchor_t> anchors = random_anchors(width, height, argsinfo.random_anchors_arg, argsinfo.random_seed_arg);
		merge_perfect_chains(anchors); // only maximal anchors!
		place_dummy_anchors(width, height, anchors);
		sort_anchors(anchors);

		// solve via ChainX precedence and plot case 2 recursions
		Image image(width, height);
		vector<long long> costs;
		vector<anchor_t> optimal_chain;
		algo::chainx_naive(anchors, algo::chaining_mode::semiglobal, costs, optimal_chain);
		plot_gap_gap_lower_diag(image, anchors, costs); // plot case 2 recursions
		plot_anchors(image, anchors); // plot all anchors
		plot_anchors(image, optimal_chain, utils::defaults::selected_anchor_color); // recolor the optimal chain
		image.writeToFile(argsinfo.debug_case_two_output_file_arg);
		assert(costs.back() == compute_chain_cost(optimal_chain, algo::chaining_mode::semiglobal));

		// solve via the would-be linearithmic solution and compare
		vector<long long> new_costs;
		algo::solve_linearithmic_debug(anchors, width, height, algo::chaining_mode::semiglobal, new_costs, costs);
		assert(new_costs.size() == costs.size() and new_costs.back() == costs.back());

		vector<anchor_t> checking_chain;
		algo::weak_backtrack(anchors, costs, algo::chaining_mode::semiglobal, checking_chain);
		cerr << "Optimal ChainX chain has cost     " << compute_chain_cost(optimal_chain, algo::chaining_mode::semiglobal) << endl;
		cerr << "Backtracked ChainX chain has cost " << compute_chain_cost(checking_chain, algo::chaining_mode::semiglobal) << endl;
		assert(compute_chain_cost(checking_chain, algo::chaining_mode::semiglobal) == compute_chain_cost(optimal_chain, algo::chaining_mode::semiglobal));
		return 0;
	}

	return 1;
}

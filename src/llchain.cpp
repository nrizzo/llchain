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

#include "command-line-parsing/cmdline.h" // cmdline_parser (gengetopt)

using std::cout, std::cerr, std::endl;
using std::string, std::to_string;
using std::vector;
using std::ifstream, std::ofstream;
using namespace llchain;
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
	if (argsinfo.chainx_flag and argsinfo.chainx_opt_flag) {
		cerr << "Error: select only one ChainX version!" << endl;
		exit(1);
	}
	if (argsinfo.chainx_original_magic_numbers_flag and not (argsinfo.chainx_flag or argsinfo.chainx_opt_flag)) {
		cerr << "Error: --chainx-original-numbers requires --chainx or --chainx-opt!" << endl;
		exit(1);
	}
	if ((argsinfo.sam_arg != NULL or argsinfo.output_arg != NULL) and argsinfo.all_to_all_flag) {
		cerr << "Error: output and all-to-all mode are incompatible!" << endl;
		exit(1);
	}

	const int anchorlength = argsinfo.anchor_length_arg;
	if (anchorlength <= 0) { cerr << "Error: pick an anchor length >= 1." << endl; exit(1); };

	algo::chaining_mode mode;
	chainx::mode chainx_mode;
	if (string(argsinfo.mode_arg) == "global") {
		mode = algo::chaining_mode::global;
		chainx_mode = chainx::mode::global;
	} else if (string(argsinfo.mode_arg) == "semiglobal") {
		mode = algo::chaining_mode::semiglobal;
		chainx_mode = chainx::mode::semiglobal;
	} else {
		cerr << "Error: pick a correct chaining mode (global/semiglobal)." << endl;
		exit(1);
	}
	cerr << "DEBUG: " << ((mode == algo::chaining_mode::global) ? "global" : "semiglobal") << " mode" << endl;

	anchor_type anchortype;
	if      (string(argsinfo.anchor_type_arg) == "MUM") anchortype = MUM;
	else if (string(argsinfo.anchor_type_arg) == "MEM") anchortype = MEM;
	else { cerr << "Error: pick a correct anchor type (MUM/MEM)." << endl; exit(1); }
	if (argsinfo.custom_anchors_arg == NULL) {
		cerr << "DEBUG: " << ((anchortype == MUM) ? "MUM" : "MEM") << " anchors of length >= " << anchorlength << endl;
	} else {
		cerr << "DEBUG: custom anchors from file " << string(argsinfo.custom_anchors_arg) << endl;
	}

	if ((argsinfo.text_arg != NULL) and (argsinfo.query_arg != NULL)) {
		vector<string> texts, text_ids; // queries, query_ids;
		kseq::read_sequences(string(argsinfo.text_arg),  texts,   text_ids);

		cerr << "DEBUG: read text sequences ";
		for (auto const &id : text_ids) cerr << id << " ";
		cerr << "of sizes ";
		for (auto const &text : texts) cerr << text.size() << " ";
		cerr << endl;

		ofstream sam_out;
		if (argsinfo.sam_arg != NULL) {
			sam_out = ofstream(string(argsinfo.sam_arg), std::ios::trunc);
			algo::write_SAM_header(sam_out);
			for (size_type t = 0; t < texts.size(); t++)
				algo::write_SAM_text(texts[t], text_ids[t], sam_out);
		}

		ofstream out;
		if (argsinfo.output_arg != NULL) {
			out = ofstream(string(argsinfo.output_arg), std::ios::trunc);
		}

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
				int chainx_revisions = -1;
				if (argsinfo.chainx_flag) {
					chainx::chainx(matches, query.size(), costs, chainx_revisions, chainx_mode, false);
				} else if (argsinfo.chainx_opt_flag) {
					chainx::chainx(matches, query.size(), costs, chainx_revisions, chainx_mode);
				} else {
					algo::weak_solve_loglinear(matches, texts[t].size(), query.size(), mode, costs);
				}
				const std::chrono::duration<double> main_chaining_time = std::chrono::steady_clock::now() - start;

				start = std::chrono::steady_clock::now();
				vector<anchor_t> chain;
				if (argsinfo.chainx_flag or argsinfo.chainx_opt_flag) {
					algo::chainx_backtrack(matches, costs, mode, chain);
				} else {
					algo::weak_backtrack(matches, costs, mode, chain);
				}
				const std::chrono::duration<double> backtrack_chaining_time = std::chrono::steady_clock::now() - start;
				const std::chrono::duration<double> query_time = std::chrono::steady_clock::now() - querystart;

				cerr << "done (" << found_anchors << " anchors, " << matches.size() - 2 << " merged, " << costs.back() << " anchored edit distance, " << seeding_time << " seeding, " << preprocessing_time << " preprocessing, " << main_chaining_time << " chaining, " << backtrack_chaining_time << " backtrack, " << query_time << " total query time" << ((argsinfo.chainx_flag or argsinfo.chainx_opt_flag) ? (", " + to_string(chainx_revisions) + " revisions") : "") << ")" << endl;

				if (argsinfo.sam_arg != NULL and chain.size() > 2) {
					algo::write_SAM_entry(texts[t], text_ids[t], query, query_id, argsinfo.store_SAM_sequence_flag, chain, mode, costs.back(), sam_out);
				}
				if (argsinfo.output_arg != NULL) {
					out << ">" << query_id << " (Reference " << text_ids[t] << ")\n";
					assert(chain.size() >= 2);
					for (size_type i = 1; i < chain.size() - 1; i++) {
						out << get<0>(chain[i])+1 << '\t' << get<1>(chain[i])+1 << '\t' << get<2>(chain[i])+1 << "\n";
					}
				}
			}
		}
		if (argsinfo.sam_arg != NULL) {
			sam_out.close();
		}
		if (argsinfo.output_arg != NULL) {
			out.close();
		}

		return 0;
	}

	if (argsinfo.all_to_all_flag) {
		vector<string> queries, query_ids;
		kseq::read_sequences(string(argsinfo.query_arg), queries, query_ids);

		vector<vector<utils::anchor_index_t>> distances(
				queries.size(),
				vector<utils::anchor_index_t>(queries.size(),
					std::numeric_limits<utils::anchor_index_t>::max())
				);
		for (size_type i = 1; i < queries.size(); i++) {
			auto start = std::chrono::steady_clock::now(), querystart = std::chrono::steady_clock::now();
			auto const index = mummer_essaMEM_wrapper::index(queries[i], anchorlength);
			const std::chrono::duration<double> index_time = std::chrono::steady_clock::now() - start;
			cerr << "DEBUG: indexed " << query_ids[i] << " in " << index_time << endl;

			for (size_type j = 0; j < i; j++) {
				cerr << "DEBUG: querying " << query_ids[j] << " in " << query_ids[i] << " (" << ((anchortype == MUM) ? "MUM" : "MEM") << " seeds of length >= " << anchorlength << ")...";
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
				int chainx_revisions = -1;
				if (argsinfo.chainx_flag) {
					chainx::chainx(matches, queries[j].size(), costs, chainx_revisions, chainx_mode, false);
				} else if (argsinfo.chainx_opt_flag) {
					chainx::chainx(matches, queries[j].size(), costs, chainx_revisions, chainx_mode);
				} else {
					algo::weak_solve_loglinear(matches, queries[i].size(), queries[j].size(), mode, costs);
				}
				const std::chrono::duration<double> main_chaining_time = std::chrono::steady_clock::now() - start;
				const std::chrono::duration<double> query_time = std::chrono::steady_clock::now() - querystart;
				distances[i][j] = costs.back();

				cerr << "done (" << found_anchors << " anchors, " << matches.size() - 2 << " merged, " << costs.back() << " anchored edit distance, " << seeding_time << " seeding, " << preprocessing_time << " preprocessing, " << main_chaining_time << " chaining, " << query_time << " total query time" << ((argsinfo.chainx_flag or argsinfo.chainx_opt_flag) ? (", " + to_string(chainx_revisions) + " revisions") : "") << ")" << endl;
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
		vector<anchor_t> chainx_chain;
		algo::chainx_solve_naive(anchors, mode, costs, chainx_chain);
		plot_gap_gap_lower_diag(image, anchors, costs); // plot case 2 recursions
		plot_anchors(image, anchors); // plot all anchors
		plot_anchors(image, chainx_chain, utils::defaults::selected_anchor_color); // recolor the optimal chain
		image.writeToFile(argsinfo.debug_case_two_output_file_arg);
		assert(costs.back() == algo::compute_chain_cost(chainx_chain, mode));

		// solve via the would-be linearithmic solution and compare
		vector<long long> new_costs;
		algo::weak_solve_loglinear_debug(anchors, width, height, mode, new_costs, costs);
		assert(new_costs.size() == costs.size() and new_costs.back() == costs.back());

		vector<anchor_t> chain;
		algo::weak_backtrack(anchors, costs, mode, chain);
		cerr << "Optimal chain has cost     " << costs.back() << endl;
		cerr << "Backtracked chain has cost " << algo::compute_chain_cost(chain, mode) << endl;
		assert(algo::compute_chain_cost(chain, mode) == algo::compute_chain_cost(chainx_chain, mode));

		// output chain
		if (argsinfo.output_arg != NULL) {
			ofstream out(string(argsinfo.output_arg), std::ios::trunc);
			out << ">" << "query" << '\n';
			assert(chain.size() >= 2);
			for (size_type i = 1; i < chain.size() - 1; i++) {
				out << get<0>(chain[i])+1 << '\t' << get<1>(chain[i])+1 << '\t' << get<2>(chain[i])+1 << "\n";
			}
			out.close();
		}

		// output SAM
		if (argsinfo.sam_arg != NULL) {
			ofstream out(string(argsinfo.sam_arg ), std::ios::trunc);
			algo::write_SAM_header(out);
			algo::write_SAM_text(string(width, 'A'), "ref", out);
			algo::write_SAM_entry(string(width, 'A'), "ref", string(height, 'A'), "query", argsinfo.store_SAM_sequence_flag, chain, mode, costs.back(), out);
			out.close();
		}
		return 0;
	}

	return 1;
}

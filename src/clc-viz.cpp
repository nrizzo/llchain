#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <limits>
#include <algorithm>
#include <random>
#include <map>
#include <cassert>
#include <set>

import Utils;
import Algo;
#include "command-line-parsing/cmdline.h" // gengetopt-generated parser

using std::cerr, std::endl, std::vector, std::string, std::tuple, std::get, std::sort, std::min, std::min_element, std::max, std::abs, std::pair, std::multiset;

int main(int argc, char **argv) 
{
	// setup
	gengetopt_args_info argsinfo;
	if (cmdline_parser(argc, argv, &argsinfo) != 0) exit(1);
	if (argsinfo.anchors_arg == 5) cerr << "Running with default option \"-n 5\"" << endl;
	if (string(argsinfo.gap_gap_lower_diagonal_output_file_arg) == "gap-gap-ld.bmp") cerr << "Running with default option \"-g gap-gap-ld.bmp\"" << endl;
	if (argsinfo.random_seed_arg == -1) cerr << "Running with default option \"-r -1\"" << endl;

	const int width = 400;
	const int height = 200;

	// anchor generation
	vector<anchor_t> anchors = random_anchors(width, height, argsinfo.anchors_arg, argsinfo.random_seed_arg);
	anchors = filter_perfect_chains(anchors); // only maximal anchors!
	place_dummy_anchors(width, height, anchors);
	sort_anchors(anchors);

	// solve via ChainX precedence and plot case 2 recursions
	Image image(width, height);
	vector<long long> costs;
	vector<anchor_t> optimal_chain;
	chainx_global_naive(anchors, costs, optimal_chain);
	plot_gap_gap_lower_diag(image, anchors, costs); // plot case 2 recursions
	plot_anchors(image, anchors); // plot all anchors
	plot_anchors(image, optimal_chain, defaults::selected_anchor_color); // recolor the optimal chain
	image.writeToFile(argsinfo.gap_gap_lower_diagonal_output_file_arg);

	// solve via the would-be linearithmic solution and compare
	vector<long long> new_costs;
	solve_global_linearithmic_naive(anchors, new_costs, optimal_chain, costs);
	assert(new_costs.size() == costs.size() and new_costs.back() == costs.back());

	return 0;
}

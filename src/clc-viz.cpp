#include <iostream>
#include <string>
#include <vector>
#include <cassert>

import utils;
import algo;
#include "command-line-parsing/cmdline.h" // gengetopt-generated parser

using std::cerr, std::endl;
using std::string;
using std::vector;
using utils::anchor_t, utils::random_anchors, utils::plot_gap_gap_lower_diag, utils::plot_anchors, utils::Image, utils::filter_perfect_chains, utils::place_dummy_anchors, utils::sort_anchors;

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
	algo::chainx_global_naive(anchors, costs, optimal_chain);
	plot_gap_gap_lower_diag(image, anchors, costs); // plot case 2 recursions
	plot_anchors(image, anchors); // plot all anchors
	plot_anchors(image, optimal_chain, utils::defaults::selected_anchor_color); // recolor the optimal chain
	image.writeToFile(argsinfo.gap_gap_lower_diagonal_output_file_arg);

	// solve via the would-be linearithmic solution and compare
	vector<long long> new_costs;
	algo::solve_global_linearithmic_naive(anchors, new_costs, optimal_chain, costs);
	assert(new_costs.size() == costs.size() and new_costs.back() == costs.back());

	return 0;
}

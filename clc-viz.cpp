#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <limits>
#include <algorithm>
#include <random>

#include "ext/grid_to_bmp.hpp"
#include "command-line-parsing/cmdline.h" // gengetopt-generated parser

using std::cerr, std::endl, std::vector, std::string, std::tuple, std::get, grid_to_bmp::Color, std::sort, std::min, std::max, std::abs;

const Color selected_anchor_color(255, 0, 0);
const Color anchor_color(0, 0, 0);
const vector<Color> palette = {
	Color(213,  94,   0), // Rusty orange
	Color(  0, 158, 115), // Bluish green
	Color(  0, 114, 178), // Bluish
	Color(240, 228,  66), // Dirty yellow
	Color(230, 159,   0), // Squash
	Color( 86, 180, 233), // Sky
	Color(204, 121, 167), // Dull pink
	Color(225, 156, 137), // Pinkish tan
	Color(137, 189, 165), // Pale teal
	Color(137, 165, 201), // Slate blue
	Color(244, 235, 147), // Dark cream
	Color(237, 189, 137), // Fawn
	Color(139,  61,   0), // Sienna
	Color(  0, 103,  75), // Spruce
	Color(  0,  74, 116), // Deep sea blue
	Color(185, 173,  14), // Muddy yellow
	Color(150, 103,   0), // Poo brown
};


typedef tuple<long long, long long, long long> anchor_t;
vector<anchor_t> random_anchors(long long width, long long height, int n, int length, int random_seed);
void place_dummy_anchors(int width, int height, vector<anchor_t> &anchors);
void plot_anchors(grid_to_bmp::BmpImage &image, vector<anchor_t> &anchors, Color color = anchor_color);
vector<tuple<long long,long long,long long>> solve_and_plot_gap_gap_forward_global(grid_to_bmp::BmpImage &image, vector<anchor_t> &anchors);

int main(int argc, char **argv) 
{
	gengetopt_args_info argsinfo;
	if (cmdline_parser(argc, argv, &argsinfo) != 0) exit(1);
	if (argsinfo.anchors_arg == 3) cerr << "Running with default option \"-n 3\"" << endl;
	if (string(argsinfo.gap_gap_output_file_arg) == "gap-gap.bmp") cerr << "Running with default option \"-g gap-gap.bmp\"" << endl;
	if (argsinfo.random_seed_arg == -1) cerr << "Running with default option \"-r -1\"" << endl;

	const int width = 200;
	const int height = 100;
	const int anchor_length = 10;

	vector<anchor_t> anchors = random_anchors(width, height, argsinfo.anchors_arg, anchor_length, argsinfo.random_seed_arg);
	place_dummy_anchors(width, height, anchors);
	std::sort (anchors.begin(), anchors.end(),
			[](const anchor_t& a,
				const anchor_t& b) -> bool
			{
			return get<0>(a) < get<0>(b);
			});

	grid_to_bmp::BmpImage image(width, height);
	vector<anchor_t> optimal_chain = solve_and_plot_gap_gap_forward_global(image, anchors);
	plot_anchors(image, anchors);
	plot_anchors(image, optimal_chain, selected_anchor_color);
	image.writeToFile(argsinfo.gap_gap_output_file_arg);
}

vector<anchor_t> random_anchors(long long width, long long height, int n, int length, int random_seed)
{
	vector<anchor_t> res;
	std::uniform_int_distribution<int> qprng(0, width  - length);
	std::uniform_int_distribution<int> tprng(0, height - length);
	std::random_device r;
	std::mt19937 gen;
	if (random_seed == -1) {
		gen.seed(r());
	} else {
		gen.seed(random_seed);
	}

	for (int k = 0; k < n; k++) {
		int qbeg = qprng(gen);
		int tbeg = tprng(gen);
		res.emplace_back(qbeg, tbeg, length);
	}


	return res;
}

void plot_anchors(grid_to_bmp::BmpImage &image, vector<anchor_t> &anchors, Color color)
{
	for (auto &a : anchors) {
		long long qbeg = get<0>(a);
		long long tbeg = get<1>(a);
		long long match_length = get<2>(a);

		if (qbeg < 0 or qbeg >= image.width or tbeg < 0 or tbeg >= image.height)
			continue;
		for (int k = 0; k < match_length; k++) {
			image.putpixel(qbeg + k, tbeg + k, color);
		}
	}
}

void place_dummy_anchors(int width, int height, vector<anchor_t> &anchors)
{
	anchors.emplace_back(-1,-1,1);
	anchors.emplace_back(width, height, 1);
}

vector<anchor_t> solve_and_plot_gap_gap_forward_global(grid_to_bmp::BmpImage &image, vector<anchor_t> &anchors)
{
	long long n = anchors.size();
	vector<long long> costs(n, 0);
	vector<long long> backtracks(n, -1);
	vector<anchor_t> optimal_chain;

	// 1. compute optimal costs (ChainX precedence)
	for(long long j = 1; j < n; j++) // ignore first dummy anchors
	{
		//compute cost[i] here
		long long find_min_cost = std::numeric_limits<long long>::max();
		long long backtrack = -1;

		long long j_a = get<0>(anchors[j]);
		long long j_b = get<0>(anchors[j]) + get<2>(anchors[j]) - 1;
		long long j_c = get<1>(anchors[j]);
		long long j_d = get<1>(anchors[j]) + get<2>(anchors[j]) - 1;

		// anchor i < anchor j 
		for(long long i = j - 1; i >= 0; i--)
		{
			long long i_a = get<0>(anchors[i]);
			long long i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
			long long i_c = get<1>(anchors[i]);
			long long i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

			if (costs[i] < std::numeric_limits<long long>::max() && i_a < j_a && i_b < j_b && i_c < j_c && i_d < j_d)
			{
				long long gap1 = max((long long)0, j_a - i_b - 1);
				long long gap2 = max((long long)0, j_c - i_d - 1);
				long long g = max(gap1,gap2);

				long long overlap1 = max((long long)0, i_b - j_a + 1);
				long long overlap2 = max((long long)0, i_d - j_c + 1);
				long long o = abs(overlap1 - overlap2);

				if (costs[i] + g + o < find_min_cost) backtrack = i;
				find_min_cost = min(find_min_cost, costs[i] + g + o);
			}
		}
		//save optimal cost at offset j
		costs[j] = find_min_cost;
		backtracks[j] = backtrack;
	}

	// 2. project and plot optimal costs
	vector<vector<long long>> gap_gap_forward_costs(
			image.width,
			vector<long long>(image.height, std::numeric_limits<long long>::max())
		);
	for(long long i = 1; i < n - 1; i++) { // ignore first and last dummy anchors
		long long i_a = get<0>(anchors[i]);
		long long i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
		long long i_c = get<1>(anchors[i]);
		long long i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

		for (long long j_a = i_b + 1; j_a < image.width; j_a++) {
			for (long long j_c = i_d + 1; j_c < image.height; j_c++) {
				long long gap1 = max((long long)0, j_a - i_b - 1);
				long long gap2 = max((long long)0, j_c - i_d - 1);
				long long g = max(gap1,gap2);
				long long gap_gap_cost = costs[i] + g;

				if (gap_gap_cost <= gap_gap_forward_costs[j_a][j_c]) { // TODO investigate ties
					gap_gap_forward_costs[j_a][j_c] = gap_gap_cost;
					image.putpixel(j_a, j_c, palette[(i + palette.size()-1) % palette.size()]);
				}
			}
		}
	}

	// 3. return optimal chain
	for (long long i = n - 1; backtracks[i] >= 0; i = backtracks[i]) {
		optimal_chain.push_back(anchors[i]);
	}
	optimal_chain.push_back(anchors[0]);
	std::reverse(optimal_chain.begin(), optimal_chain.end());
	return optimal_chain;
}

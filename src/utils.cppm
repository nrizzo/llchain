module; // global module fragment

#include <vector>
#include <tuple>
#include <random>
#include <algorithm> // std::sort
#include <iostream>
#include "grid_to_bmp.hpp"

export module utils; // this module

using std::vector;
using std::tuple;
using std::uniform_int_distribution, std::random_device, std::mt19937; // random
using std::max; // algorithm
using std::cerr, std::endl;
using grid_to_bmp::Color, grid_to_bmp::BmpImage; // grid_to_bmp

namespace utils {
export typedef grid_to_bmp::BmpImage Image;

export namespace defaults
{
	// basic colors and a palette
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

	// constants for anchor generation
	const int anchor_min_length = 1;
	const int anchor_max_length = 30;
} // namespace defaults

/* 
 * an anchor is represented by triple
 * (0-based start in T, 0-based start in Q, length)
 * NB: indices can be negative/greater than T-1, Q-1 (see place_dummy_anchors)
 */
export typedef long long anchor_index_t;
export typedef tuple<anchor_index_t, anchor_index_t, anchor_index_t> anchor_t;

/*
 * computes cost of connecting anc2 after anc1 (max gap + absolute overlap diff)
 */
export anchor_index_t connect(const anchor_t &anc1, const anchor_t &anc2);
export inline anchor_index_t connect(anchor_index_t a1, anchor_index_t b1, anchor_index_t c1, anchor_index_t d1, anchor_index_t a2, anchor_index_t b2, anchor_index_t c2, anchor_index_t d2);

/*
 * generate n random anchors
 * NB: anchors can overlap, intersect, be contained in one another
 * (see filter_perfect_chains)
 */
export vector<anchor_t> random_anchors(
		anchor_index_t Tlength,
		anchor_index_t Qlength,
		int n,
		int random_seed
	);

/*
 * return a subset of the input anchors that does not contain perfect chains of
 * length >1 (connect(a_i, a_j) == 0)
 */
export vector<anchor_t> filter_perfect_chains(const vector<anchor_t> &anchors);

/*
 * adds dummy anchors a_start and a_end to the end of anchors
 * (see sort_anchors)
 */
export void place_dummy_anchors(
		anchor_index_t Tlength,
		int Qlength,
		vector<anchor_t> &anchors);

/*
 * sorts anchors according to starting position in T
 */
export void sort_anchors(vector<anchor_t> &anchors);

/*
 * draws the given anchors on top of the image
 */
export void plot_anchors(
		BmpImage &image,
		const vector<anchor_t> &anchors,
		Color color = defaults::anchor_color
	);

/*
 * plots the optimal case-2 recursions forward
 * NB: assumes anchors/costs to be sorted by sort_anchors
 * NB: takes O(Tlength x Qlength) space and time
 * NB: ties are broken by picking most recent anchor
 */
export void plot_gap_gap_lower_diag(
		BmpImage &image,
		vector<anchor_t> &anchors,
		vector<anchor_index_t> &costs
	);

export
inline
anchor_index_t connect(const anchor_t &anc1, const anchor_t &anc2)
{
	const anchor_index_t a1 = get<0>(anc1);
	const anchor_index_t b1 = get<0>(anc1) + get<2>(anc1);
	const anchor_index_t c1 = get<1>(anc1);
	const anchor_index_t d1 = get<1>(anc1) + get<2>(anc1);

	const anchor_index_t a2 = get<0>(anc2);
	const anchor_index_t b2 = get<0>(anc2) + get<2>(anc2);
	const anchor_index_t c2 = get<1>(anc2);
	const anchor_index_t d2 = get<1>(anc2) + get<2>(anc2);

	const anchor_index_t Tgap = max((anchor_index_t)0, a2 - b1);
	const anchor_index_t Qgap = max((anchor_index_t)0, c2 - d1);
	const anchor_index_t g = max(Tgap, Qgap);
	const anchor_index_t Tovl = max((anchor_index_t)0, b1 - a2);
	const anchor_index_t Qovl = max((anchor_index_t)0, d1 - c2);
	const anchor_index_t o = abs(Tovl - Qovl);

	return g + o;
}

export
inline
anchor_index_t connect(anchor_index_t a1, anchor_index_t b1, anchor_index_t c1, anchor_index_t d1, anchor_index_t a2, anchor_index_t b2, anchor_index_t c2, anchor_index_t d2)
{
	const anchor_index_t Tgap = max((anchor_index_t)0, a2 - b1);
	const anchor_index_t Qgap = max((anchor_index_t)0, c2 - d1);
	const anchor_index_t g = max(Tgap, Qgap);
	const anchor_index_t Tovl = max((anchor_index_t)0, b1 - a2);
	const anchor_index_t Qovl = max((anchor_index_t)0, d1 - c2);
	const anchor_index_t o = abs(Tovl - Qovl);

	return g + o;
}

export vector<anchor_t> random_anchors(anchor_index_t width, anchor_index_t height, int n, int random_seed)
{
	vector<anchor_t> res;
	uniform_int_distribution<int> qprng(0, (width  - defaults::anchor_max_length) * 75 / 100);
	uniform_int_distribution<int> tprng(0, (height - defaults::anchor_max_length) * 75 / 100);
	uniform_int_distribution<int> lprng(1, defaults::anchor_max_length - defaults::anchor_min_length);
	random_device r;
	mt19937 gen;

	if (random_seed == -1) {
		int random_seed = r();
		cerr << "Running with random seed " << random_seed << endl;
		gen.seed(random_seed);
	} else {
		gen.seed(random_seed);
	}

	for (int k = 0; k < n; k++) {
		int qbeg = qprng(gen);
		int tbeg = tprng(gen);
		int length = lprng(gen);
		res.emplace_back(qbeg, tbeg, length);
	}

	return res;
}

export vector<anchor_t> filter_perfect_chains(vector<anchor_t> &anchors)
{
	// sort anchors by diagonal first, start position second
	std::sort(anchors.begin(), anchors.end(),
			[](const anchor_t& a,
				const anchor_t& b) -> bool
			{
				anchor_index_t diag_a = get<0>(a) - get<1>(a);
				anchor_index_t diag_b = get<0>(b) - get<1>(b);
				return (diag_a < diag_b or
						(diag_a == diag_b and
						 get<0>(a) < get<0>(b)));
			});

	vector<anchor_t> solution;
	for (anchor_index_t i = 0; i < anchors.size(); i++) {
		if (solution.size() > 0) {
			anchor_index_t i_diag = get<0>(anchors[i]) - get<1>(anchors[i]);
			anchor_index_t i_a = get<0>(anchors[i]);
			anchor_index_t s_diag = get<0>(solution.back()) - get<1>(solution.back());
			anchor_index_t s_b = get<0>(solution.back()) + get<2>(solution.back());

			if (i_diag == s_diag and i_a <= s_b) {
				// do not add, connect(i,s) == 0
			} else {
				solution.push_back(anchors[i]);
			}
		} else {
			solution.push_back(anchors[i]);
		}
	}
	return solution;
}

export void sort_anchors(vector<anchor_t> &anchors)
{
	std::sort(anchors.begin(), anchors.end(),
			[](const anchor_t& a,
				const anchor_t& b) -> bool
			{
			return get<0>(a) < get<0>(b);
			});
}

export void plot_anchors(BmpImage &image, const vector<anchor_t> &anchors, Color color)
{
	for (const auto &a : anchors) {
		anchor_index_t qbeg = get<0>(a);
		anchor_index_t tbeg = get<1>(a);
		anchor_index_t match_length = get<2>(a);

		if (qbeg < 0 or qbeg + match_length > image.width or
		    tbeg < 0 or tbeg + match_length > image.height  )
			continue;
		for (int k = 0; k < match_length; k++) {
			image.putpixel(qbeg + k, tbeg + k, color);
		}
	}
}

export void place_dummy_anchors(int width, int height, vector<anchor_t> &anchors)
{
	anchors.emplace_back(-1,-1,1);
	anchors.emplace_back(width, height, 1);
}

export void plot_gap_gap_lower_diag(BmpImage &image, vector<anchor_t> &anchors, vector<anchor_index_t> &costs)
{
	anchor_index_t n = anchors.size();

	// project and plot optimal costs
	vector<vector<anchor_index_t>> gap_gap_forward_costs(
			image.width,
			vector<anchor_index_t>(image.height, std::numeric_limits<anchor_index_t>::max())
		);
	for(anchor_index_t i = 1; i < n - 1; i++) { // ignore first and last dummy anchors
		anchor_index_t i_a = get<0>(anchors[i]);
		anchor_index_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
		anchor_index_t i_c = get<1>(anchors[i]);
		anchor_index_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

		for (anchor_index_t j_a = i_b + 1; j_a < image.width; j_a++) {
			for (anchor_index_t j_c = i_d + 1; j_c < image.height; j_c++) {
				anchor_index_t gap1 = max((anchor_index_t)0, j_a - i_b - 1);
				anchor_index_t gap2 = max((anchor_index_t)0, j_c - i_d - 1);
				anchor_index_t g = max(gap1,gap2);
				anchor_index_t gap_gap_cost = costs[i] + g;

				if (i_a - i_c <= j_a - j_c and
						gap_gap_cost <= gap_gap_forward_costs[j_a][j_c]) {
					gap_gap_forward_costs[j_a][j_c] = gap_gap_cost;
					image.putpixel(
							j_a,
							j_c,
							defaults::palette[
								(i + defaults::palette.size()-1) %
								defaults::palette.size()]
					      );
				}
			}
		}
	}
}
} // namespace utils

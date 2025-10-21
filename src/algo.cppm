module; // global module fragment

#include <vector>
#include <tuple>
#include <algorithm> // std::sort
#include <iostream>
#include <map>
#include <set>
#include <cassert>

export module algo; // this module

import utils;
import MinSegmentTree;

using std::vector;
using std::tuple, std::pair;
using std::max, std::min, std::abs;
using std::cerr, std::endl;
using std::map;
using std::multiset;
using utils::anchor_t, utils::connect;
typedef utils::anchor_index_t ai_t;

namespace algo {
/*
 * solves chainx-precedence colinear chaining via DP, returns an optimal chain
 *   via backtracking; adapted from github.com/at-cg/ChainX
 * NB: O(n^2) time
 * NB: assumes anchors contains dummies (see place_dummy_anchors)
 * NB: assumes sorted anchors (see sort_anchors)
 */
export void chainx_global_naive(
		const vector<anchor_t> &anchors,
		vector<ai_t> &costs_out,
		vector<anchor_t> &chain_out
	);

/*
 * solves (relaxed chainx-precedence) colinear chaining via DP and mimicks the
 *   solution/data structures of our linearithmic-time prototype
 * NB: O(n^2) time or worse
 * NB: assumes anchors contains dummies (see place_dummy_anchors)
 * NB: assumes sorted anchors (see sort_anchors)
 */
export void solve_global_linearithmic_naive(
		const vector<anchor_t> &anchors,
		const ai_t Tlength,
		const ai_t Qlength,
		vector<ai_t> &costs_out,
		vector<anchor_t> &chain_out,
		const vector<ai_t> &correct_costs
	);

// start of implementation
export
void chainx_global_naive(
		const vector<anchor_t> &anchors,
		vector<ai_t> &costs_out,
		vector<anchor_t> &chain_out
) {
	const ai_t n = anchors.size();
	costs_out = vector<ai_t>(n); // partial chainx-prec DP costs
	vector<ai_t> backtracks(n, -1); // index of an optimal previous anchor
	chain_out.clear();

	costs_out[0] = 0; // first dummy anchor
	for (ai_t i = 1; i < n; i++) {
		ai_t i_cost = std::numeric_limits<ai_t>::max();
		ai_t backtrack = -1;

		const ai_t i_a = get<0>(anchors[i]);
		const ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
		const ai_t i_c = get<1>(anchors[i]);
		const ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

		for (ai_t j = i - 1; j >= 0; j--) { // anchor j < anchor i 
			const ai_t j_a = get<0>(anchors[j]);
			const ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]);
			const ai_t j_c = get<1>(anchors[j]);
			const ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]);

			if (costs_out[j] < std::numeric_limits<ai_t>::max() and
			    j_a < i_a and j_b < i_b and j_c < i_c and j_d < i_d) { // chainx precedence
				const ai_t c = costs_out[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d);
				if (c < i_cost) backtrack = j;
				i_cost = min(i_cost, c);
			}
		}
		//save optimal cost at offset i
		costs_out[i] = i_cost;
		backtracks[i] = backtrack;
	}

	// trace back an optimal chain
	for (ai_t i = n - 1; backtracks[i] >= 0; i = backtracks[i]) {
		chain_out.push_back(anchors[i]);
	}
	chain_out.push_back(anchors[0]);
	std::reverse(chain_out.begin(), chain_out.end());
}

/*
 * case one: gap-gap case, bigger gap in Q
 *   requires a 1-dimensional range minimum query data structure where the
 *   (static) positions are anchor diagonals and the values (semi-dynamic,
 *   better values only) can be negative
 */
struct case_one_index {
	MinSegmentTree<ai_t> recursive_values; // diagonal -> min(C[j] - q_e^j)
};
struct case_one_index init_case_one(ai_t Tlength, ai_t Qlength);
ai_t compute_case_one(const case_one_index &I, const anchor_t &a_i);
void update_startpoint_case_one(case_one_index &I, const anchor_t &a_j, const ai_t j_cost);
void  update_endpoint_case_one (case_one_index &I, const anchor_t &a_j, const ai_t j_cost);

/*
 * help function that checks all recursive cases considered by case one
 */
ai_t compute_case_one_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs);

// unused functions below
struct case_one_index_naive {
	std::map<ai_t, ai_t> recursive_values; // diagonal -> min(C[j] - q_e^j)
};
struct case_one_index_naive init_case_one_naive();
ai_t compute_case_one_naive(const case_one_index_naive &I, const anchor_t &a_i);
void update_startpoint_case_one_naive(case_one_index_naive &I, const anchor_t &a_j, const ai_t j_cost);
void  update_endpoint_case_one_naive (case_one_index_naive &I, const anchor_t &a_j, const ai_t j_cost);

enum line
{
	horizontal,
	diagonal
};

struct case_two_index {
	// TODO implement an efficient version
	const vector<anchor_t> &anchors;
	vector<pair<ai_t,line>> delimiting_lines; // element is pair of anchor index and line type
	vector<ai_t> optimal_recursion_values;
};

struct case_two_index init_case_two(const vector<anchor_t> &anchors) {
	case_two_index I = {anchors, {{0, horizontal}, {anchors.size()-1,horizontal}}, {std::numeric_limits<ai_t>::max()}};
	return I;
}

ai_t get_c(case_two_index &I, ai_t line_index, ai_t current_a)
{
	assert(line_index >= 0 and line_index < I.delimiting_lines.size());
	ai_t i = get<0>(I.delimiting_lines[line_index]);
	if (get<1>(I.delimiting_lines[line_index]) == horizontal) {
		ai_t i_d = get<1>(I.anchors[i]) + get<2>(I.anchors[i]);
		return i_d;
	} else { // diagonal
		ai_t i_diag = get<0>(I.anchors[i]) - get<1>(I.anchors[i]);
		return current_a - i_diag;
	}
}

void prune_shadowed_delimiting_lines(case_two_index &I, ai_t sweeping_line_a)
{
	if (I.delimiting_lines.size() > 0) {
		assert(I.delimiting_lines.size() == I.optimal_recursion_values.size() + 1);
		assert(I.delimiting_lines[0].second == horizontal);
		for(ai_t l = 1; l < I.delimiting_lines.size() - 1; l++) { // l>=2?
			// has a line reached the following one?
			pair<ai_t,line> l1 = I.delimiting_lines[l];
			pair<ai_t,line> l2 = I.delimiting_lines[l + 1];
			if (l1.second == diagonal and l2.second == horizontal) {
				ai_t l1_diag = get<0>(I.anchors[l1.first]) - get<1>(I.anchors[l1.first]);
				ai_t l2_d = get<1>(I.anchors[l2.first]) + get<2>(I.anchors[l2.first]);

				if (sweeping_line_a - l1_diag >= l2_d) {
					if (l == I.delimiting_lines.size() - 2 or I.optimal_recursion_values[l-1] <= I.optimal_recursion_values[l+1]) {
						// diagonal shadows horizontal
						cerr << "DEBUG: line " << I.delimiting_lines[l].first << ((I.delimiting_lines[l].second == horizontal) ? "hor" : "diag") << " shadows " << I.delimiting_lines[l+1].first << ((I.delimiting_lines[l+1].second == horizontal) ? "hor" : "diag") << endl;
						I.delimiting_lines.erase(I.delimiting_lines.begin() + l + 1);
						I.optimal_recursion_values.erase(I.optimal_recursion_values.begin() + l);
						l -= 1;
					} else {
						// horizontal shadows diagonal
						cerr << "DEBUG: line " << I.delimiting_lines[l].first << ((I.delimiting_lines[l].second == horizontal) ? "hor" : "diag") << " is shadowed by " << I.delimiting_lines[l+1].first << ((I.delimiting_lines[l+1].second == horizontal) ? "hor" : "diag") << endl;
						I.delimiting_lines.erase(I.delimiting_lines.begin() + l);
						I.optimal_recursion_values.erase(I.optimal_recursion_values.begin() + l);
						l = max(l - 2, 0LL);
					}
				}
			}
		}
	}
}

ai_t compute_case_two(case_two_index &I, const anchor_t &anchor_i)
{
	// TODO binary search
	ai_t i_a = get<0>(anchor_i);
	ai_t i_c = get<1>(anchor_i);
	prune_shadowed_delimiting_lines(I, i_a);

	// no recursion available
	if (I.delimiting_lines.size() == 0 or
			i_c < get_c(I, 0, i_a) or
			i_c > get_c(I, I.delimiting_lines.size() - 1, i_a))
		return std::numeric_limits<ai_t>::max();

	// recursion available
	ai_t i = 0;
	while (i_c >= get_c(I, i, i_a)) i++;
	if (I.delimiting_lines[i-1].second == diagonal and i_c == get_c(I, i-1, i_a)) i--;
	if (I.optimal_recursion_values[i-1] != std::numeric_limits<ai_t>::max()) {
		return I.optimal_recursion_values[i-1] + i_a;
	} else {
		return std::numeric_limits<ai_t>::max();
	}
}

void update_startpoint_case_two(case_two_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	// do nothing
}

void update_endpoint_case_two(case_two_index &I, ai_t j, const anchor_t &a_j, const ai_t j_cost)
{
	ai_t j_b = get<0>(a_j) + get<2>(a_j);
	ai_t j_d = get<1>(a_j) + get<2>(a_j);
	prune_shadowed_delimiting_lines(I, j_b);

	// value of the two lines to potentially insert
	ai_t recursion_value = j_cost - j_b;

	// find position in array
	ai_t i = 0;
	while (i < I.delimiting_lines.size() and j_d > get_c(I, i, j_b)) i++;

	if (I.delimiting_lines.size() == 0) {
		I.optimal_recursion_values.insert(
				I.optimal_recursion_values.begin(),
				recursion_value);
		I.delimiting_lines.insert(I.delimiting_lines.begin(), std::make_pair(j, diagonal));
		I.delimiting_lines.insert(I.delimiting_lines.begin(), std::make_pair(j, horizontal));
	} else if (i == 0 or i == I.delimiting_lines.size()) {
		I.optimal_recursion_values.insert(
				I.optimal_recursion_values.begin() + i,
				recursion_value);
		I.optimal_recursion_values.insert(
				I.optimal_recursion_values.begin() + i,
				recursion_value);
		I.delimiting_lines.insert(I.delimiting_lines.begin() + i, std::make_pair(j, diagonal));
		I.delimiting_lines.insert(I.delimiting_lines.begin() + i, std::make_pair(j, horizontal));
	} else {
		if (j_d < get_c(I, i, j_b)) {
			// case 1: no line intersects with the two new ones
			if (recursion_value <= I.optimal_recursion_values[i-1]) {
				I.optimal_recursion_values.insert(
						I.optimal_recursion_values.begin() + i,
						I.optimal_recursion_values[i-1]);
				I.optimal_recursion_values.insert(
						I.optimal_recursion_values.begin() + i,
						recursion_value);
				I.delimiting_lines.insert(I.delimiting_lines.begin() + i, std::make_pair(j, diagonal));
				I.delimiting_lines.insert(I.delimiting_lines.begin() + i, std::make_pair(j, horizontal));
			}
			// else this value is not optimal and is immediately shadowed
		} else {
			// case 2: some line intersects with the two new ones
			assert(j_d == get_c(I, i, j_b)); // sanity check
			assert(i >= I.delimiting_lines.size() - 1 or j_d != get_c(I, i+1, j_b)); // two adjacent lines should never have the same coordinate if the anchors are maximal
			if (j_d == get_c(I, i, j_b) and recursion_value <= I.optimal_recursion_values[i]) {
				if (get<1>(I.delimiting_lines[i]) == horizontal) {
					I.delimiting_lines.erase(I.delimiting_lines.begin() + i);
					I.delimiting_lines.insert(I.delimiting_lines.begin() + i, std::make_pair(j, diagonal));
					I.delimiting_lines.insert(I.delimiting_lines.begin() + i, std::make_pair(j, horizontal));
					I.optimal_recursion_values.insert(
							I.optimal_recursion_values.begin() + i,
							recursion_value);
				} else if (get<1>(I.delimiting_lines[i]) == diagonal) {
					I.delimiting_lines.erase(I.delimiting_lines.begin() + i);
					I.delimiting_lines.insert(I.delimiting_lines.begin() + i, std::make_pair(j, diagonal));
					I.delimiting_lines.insert(I.delimiting_lines.begin() + i, std::make_pair(j, horizontal));
					I.optimal_recursion_values.insert(
							I.optimal_recursion_values.begin() + i,
							recursion_value);
				} else {
					assert(false);
				}
			}
		}
	}
}

ai_t naive_compute_case_two(const vector<anchor_t> &anchors, ai_t i, const vector<ai_t> &costs_out)
{
	// naively compute what we THINK we are computing in case two
	ai_t find_min_cost = std::numeric_limits<ai_t>::max();

	ai_t i_a = get<0>(anchors[i]);
	ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
	ai_t i_c = get<1>(anchors[i]);
	ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

	// anchor j < anchor i 
	for(ai_t j = i - 1; j > 0; j--) // dummy anchor[0] is treated separately!
	{
		ai_t j_a = get<0>(anchors[j]);
		ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]) - 1;
		ai_t j_c = get<1>(anchors[j]);
		ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]) - 1;

		if (costs_out[j] < std::numeric_limits<ai_t>::max() and 
				j_b < i_a and j_d < i_c and j_a - j_c <= i_a - i_c) // case 2
		{
			ai_t gap1 = max((ai_t)0, i_a - j_b - 1);
			ai_t gap2 = max((ai_t)0, i_c - j_d - 1);
			ai_t g = max(gap1,gap2);

			ai_t overlap1 = max((ai_t)0, j_b - i_a + 1);
			ai_t overlap2 = max((ai_t)0, j_d - i_c + 1);
			ai_t o = abs(overlap1 - overlap2);

			find_min_cost = min(find_min_cost, costs_out[j] + g + o);
		}
	}

	return find_min_cost;
}

struct case_three_index {
	// TODO implement an efficient version
	std::map<ai_t, multiset<ai_t>> recursive_values; // diagonal -> multiset of values C[j] + diag_j
};

struct case_three_index init_case_three()
{
	return case_three_index();
}

ai_t compute_case_three(const case_three_index &I, const anchor_t &a_i)
{
	ai_t i_a = get<0>(a_i);
	ai_t i_c = get<1>(a_i);
	ai_t i_diag = i_a - i_c;

	// naive range min query
	ai_t rec_min = std::numeric_limits<ai_t>::max();
	for (auto it = I.recursive_values.upper_bound(i_diag); it != I.recursive_values.end(); it++) {
		const multiset<ai_t> rec_values = it->second;
		rec_min = min(*min_element(rec_values.begin(), rec_values.end()), rec_min);
	}

	if (rec_min < std::numeric_limits<ai_t>::max()) {
		return rec_min - i_diag;
	} else {
		return std::numeric_limits<ai_t>::max();
	}
}

void update_startpoint_case_three(case_three_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	ai_t j_a = get<0>(a_j);
	ai_t j_c = get<1>(a_j);
	ai_t j_diag = j_a - j_c;
	ai_t rec_value = j_cost + j_diag;

	if (I.recursive_values.contains(j_diag)) {
		I.recursive_values[j_diag].insert(rec_value);
	} else {
		I.recursive_values[j_diag] = { rec_value };
	}
}

void update_endpoint_case_three(case_three_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	ai_t j_a = get<0>(a_j);
	ai_t j_c = get<1>(a_j);
	ai_t j_diag = j_a - j_c;
	ai_t rec_value = j_cost + j_diag;

	I.recursive_values[j_diag].erase(
			I.recursive_values[j_diag].find(rec_value)
		);
	if (I.recursive_values[j_diag].size() == 0) {
		I.recursive_values.erase(j_diag);
	}
}

ai_t naive_compute_case_three(const vector<anchor_t> &anchors, ai_t i, const vector<ai_t> &costs_out)
{
	// naively compute what we THINK we are computing in case three
	ai_t find_min_cost = std::numeric_limits<ai_t>::max();

	ai_t i_a = get<0>(anchors[i]);
	ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
	ai_t i_c = get<1>(anchors[i]);
	ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

	// anchor j < anchor i 
	for(ai_t j = i - 1; j > 0; j--) // dummy anchor[0] is treated separately!
	{
		ai_t j_a = get<0>(anchors[j]);
		ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]) - 1;
		ai_t j_c = get<1>(anchors[j]);
		ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]) - 1;

		if (costs_out[i] < std::numeric_limits<ai_t>::max() and 
				j_a <= i_a and i_a <= j_b and j_a - j_c > i_a - i_c) // case 3
		{
			ai_t gap1 = max((ai_t)0, i_a - j_b - 1);
			ai_t gap2 = max((ai_t)0, i_c - j_d - 1);
			ai_t g = max(gap1,gap2);

			ai_t overlap1 = max((ai_t)0, j_b - i_a + 1);
			ai_t overlap2 = max((ai_t)0, j_d - i_c + 1);
			ai_t o = abs(overlap1 - overlap2);

			find_min_cost = min(find_min_cost, costs_out[j] + g + o);
		}
	}

	return find_min_cost;
}

struct case_four_index {
	// TODO implement an efficient version
	std::map<ai_t, multiset<ai_t>> recursive_values; // diagonal -> multiset of values C[j] + diag_j
};

struct case_four_index init_case_four()
{
	return case_four_index();
}

ai_t compute_case_four(const case_four_index &I, const anchor_t &a_i)
{
	ai_t i_a = get<0>(a_i);
	ai_t i_c = get<1>(a_i);
	ai_t i_diag = i_a - i_c;

	// naive range min query
	ai_t rec_min = std::numeric_limits<ai_t>::max();
	for (auto it = I.recursive_values.begin(); it != I.recursive_values.upper_bound(i_diag); it++) {
		const multiset<ai_t> rec_values = it->second;
		rec_min = min(*min_element(rec_values.begin(), rec_values.end()), rec_min);
	}

	if (rec_min < std::numeric_limits<ai_t>::max()) {
		return rec_min + i_diag;
	} else {
		return std::numeric_limits<ai_t>::max();
	}
}

void update_startpoint_case_four(case_four_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	ai_t j_a = get<0>(a_j);
	ai_t j_c = get<1>(a_j);
	ai_t j_diag = j_a - j_c;
	ai_t rec_value = j_cost - j_diag;

	if (I.recursive_values.contains(j_diag)) {
		I.recursive_values[j_diag].insert(rec_value);
	} else {
		I.recursive_values[j_diag] = { rec_value };
	}
}

void update_endpoint_case_four(case_four_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	ai_t j_a = get<0>(a_j);
	ai_t j_c = get<1>(a_j);
	ai_t j_diag = j_a - j_c;
	ai_t rec_value = j_cost - j_diag;

	I.recursive_values[j_diag].erase(
			I.recursive_values[j_diag].find(rec_value)
		);
	if (I.recursive_values[j_diag].size() == 0) {
		I.recursive_values.erase(j_diag);
	}
}

ai_t naive_compute_case_four(const vector<anchor_t> &anchors, ai_t i, const vector<ai_t> &costs_out)
{
	// naively compute what we THINK we are computing in case four
	ai_t find_min_cost = std::numeric_limits<ai_t>::max();

	ai_t i_a = get<0>(anchors[i]);
	ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
	ai_t i_c = get<1>(anchors[i]);
	ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

	// anchor j < anchor i 
	for(ai_t j = i - 1; j > 0; j--) // dummy anchor[0] is treated separately!
	{
		ai_t j_a = get<0>(anchors[j]);
		ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]) - 1;
		ai_t j_c = get<1>(anchors[j]);
		ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]) - 1;

		if (costs_out[i] < std::numeric_limits<ai_t>::max() and 
				j_a <= i_a and i_a <= j_b and j_a - j_c <= i_a - i_c) // case 4
		{
			ai_t gap1 = max((ai_t)0, i_a - j_b - 1);
			ai_t gap2 = max((ai_t)0, i_c - j_d - 1);
			ai_t g = max(gap1,gap2);

			ai_t overlap1 = max((ai_t)0, j_b - i_a + 1);
			ai_t overlap2 = max((ai_t)0, j_d - i_c + 1);
			ai_t o = abs(overlap1 - overlap2);

			find_min_cost = min(find_min_cost, costs_out[j] + g + o);
		}
	}

	return find_min_cost;
}

struct case_five_index {
	// TODO implement an efficient version
	std::map<ai_t, ai_t> recursive_values; // j_d -> min(C[j] - diag_j)
};

struct case_five_index init_case_five()
{
	return case_five_index();
}

ai_t compute_case_five(const case_five_index &I, const anchor_t &a_i) {
	ai_t i_a = get<0>(a_i);
	ai_t i_c = get<1>(a_i);
	ai_t i_d = get<1>(a_i) + get<2>(a_i);
	ai_t i_diag = i_a - i_c;

	// naive range min query
	ai_t rec_min = std::numeric_limits<ai_t>::max();
	for (auto it = I.recursive_values.lower_bound(i_c); it != I.recursive_values.upper_bound(i_d - 1); it++) {
		const ai_t rec_value = it->second;
		rec_min = std::min(rec_value, rec_min);
	}

	if (rec_min < std::numeric_limits<ai_t>::max()) {
		return rec_min + i_diag;
	} else {
		return std::numeric_limits<ai_t>::max();
	}
}

void update_startpoint_case_five(case_five_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	// do nothing
}

void update_endpoint_case_five(case_five_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	ai_t j_a = get<0>(a_j);
	ai_t j_c = get<1>(a_j);
	ai_t j_d = get<1>(a_j) + get<2>(a_j);
	ai_t j_diag = j_a - j_c;
	ai_t rec_value = j_cost - j_diag;

	if (I.recursive_values.contains(j_d-1)) {
		I.recursive_values[j_d-1] = std::min(I.recursive_values[j_d-1], rec_value);
	} else {
		I.recursive_values[j_d-1] = rec_value;
	}
}

ai_t naive_compute_case_five(const vector<anchor_t> &anchors, ai_t i, const vector<ai_t> &costs_out)
{
	//compute cost[i] here
	ai_t find_min_cost = std::numeric_limits<ai_t>::max();

	ai_t i_a = get<0>(anchors[i]);
	ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
	ai_t i_c = get<1>(anchors[i]);
	ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

	// anchor j < anchor i 
	for(ai_t j = i - 1; j > 0; j--) // dummy anchor[0] is treated separately!
	{
		ai_t j_a = get<0>(anchors[j]);
		ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]) - 1;
		ai_t j_c = get<1>(anchors[j]);
		ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]) - 1;

		if (costs_out[i] < std::numeric_limits<ai_t>::max() and 
				//j_a < i_a and j_b < i_b and j_c < i_c and j_d < i_d and // TODO remove this?
				j_b < i_a and i_c <= j_d and j_d <= i_d) // case 5
		{
			ai_t gap1 = max((ai_t)0, i_a - j_b - 1);
			ai_t gap2 = max((ai_t)0, i_c - j_d - 1);
			ai_t g = max(gap1,gap2);

			ai_t overlap1 = max((ai_t)0, j_b - i_a + 1);
			ai_t overlap2 = max((ai_t)0, j_d - i_c + 1);
			ai_t o = abs(overlap1 - overlap2);

			find_min_cost = min(find_min_cost, costs_out[j] + g + o);
		}
	}

	return find_min_cost;
}

export void solve_global_linearithmic_naive(
		const vector<anchor_t> &anchors,
		const ai_t Tlength,
		const ai_t Qlength,
		vector<ai_t> &costs_out,
		vector<anchor_t> &chain_out,
		const vector<ai_t> &correct_costs
	)
{
	// NB anchors down here are closed-open intervals [i_a..i_b), [i_c..i_d) MOST OF THE TIME TODO
	ai_t n = anchors.size();
	costs_out = vector<ai_t>(n, 0);
	chain_out.clear();

	vector<ai_t> points; // horizontal line sweep (+j means start, -j means end)
	points.reserve(2 * n - 2);
	for (ai_t i = 1; i < n; i++) // skip starting dummy anchor
		points.push_back(-i);
	for (ai_t i = 1; i < n; i++) // skip starting dummy anchor
		points.push_back(i);
	std::stable_sort(points.begin(), points.end(),
			[&](const ai_t i,
				const ai_t j) -> bool
			{
			return (((i >= 0) ? std::get<0>(anchors[i]) : std::get<0>(anchors[-i]) + std::get<2>(anchors[-i])) <
					((j >= 0) ? std::get<0>(anchors[j]) : std::get<0>(anchors[-j]) + std::get<2>(anchors[-j])));
			});

	case_one_index   I_one   = init_case_one(Tlength, Qlength);
	case_two_index   I_two   = init_case_two(anchors);
	case_three_index I_three = init_case_three();
	case_four_index  I_four  = init_case_four();
	case_five_index  I_five  = init_case_five();

	// compute optimal costs (ChainX RELAXED precedence)
	for (ai_t point : points) {
		if (point >= 0) { // startpoint
			ai_t i = point;
			ai_t i_a = get<0>(anchors[i]);
			ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
			ai_t i_c = get<1>(anchors[i]);
			ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

			// compute cost[i]
			ai_t cost = costs_out[0] + max(i_a, i_c); // connect(a_0, a_i)
			assert(compute_case_one  (I_one,   anchors[i]) == compute_case_one_debug (anchors, i, costs_out));
			assert(compute_case_two  (I_two,   anchors[i]) == naive_compute_case_two  (anchors, i, costs_out));
			assert(compute_case_three(I_three, anchors[i]) == naive_compute_case_three(anchors, i, costs_out));
			assert(compute_case_four (I_four,  anchors[i]) == naive_compute_case_four (anchors, i, costs_out));
			assert(compute_case_five (I_five,  anchors[i]) == naive_compute_case_five (anchors, i, costs_out));

			cost = std::min(cost, compute_case_one  (I_one,   anchors[i]));
			cost = std::min(cost, compute_case_two  (I_two,   anchors[i]));
			cost = std::min(cost, compute_case_three(I_three, anchors[i]));
			cost = std::min(cost, compute_case_four (I_four,  anchors[i]));
			cost = std::min(cost, compute_case_five (I_five,  anchors[i]));
			costs_out[i] = cost;
			assert(costs_out[i] <= correct_costs[i]);

			update_startpoint_case_one  (I_one,   anchors[i], costs_out[i]);
			update_startpoint_case_two  (I_two,   anchors[i], costs_out[i]);
			update_startpoint_case_three(I_three, anchors[i], costs_out[i]);
			update_startpoint_case_four (I_four,  anchors[i], costs_out[i]);
			update_startpoint_case_five (I_five,  anchors[i], costs_out[i]);
		} else { // endpoint
			ai_t i = -point;
			update_endpoint_case_one  (I_one,    anchors[i], costs_out[i]);
			update_endpoint_case_two  (I_two, i, anchors[i], costs_out[i]);
			update_endpoint_case_three(I_three,  anchors[i], costs_out[i]);
			update_endpoint_case_four (I_four,   anchors[i], costs_out[i]);
			update_endpoint_case_five (I_five,   anchors[i], costs_out[i]);
		}
	}
}

struct case_one_index_naive init_case_one_naive()
{
	return case_one_index_naive();
}
ai_t compute_case_one_naive(const case_one_index_naive &I, const anchor_t &a_i)
{
	ai_t i_a = get<0>(a_i);
	ai_t i_c = get<1>(a_i);
	ai_t i_diag = i_a - i_c;

	// naive range min query
	ai_t rec_min = std::numeric_limits<ai_t>::max();
	for (auto it = I.recursive_values.upper_bound(i_diag); it != I.recursive_values.end(); it++) {
		const ai_t rec_value = it->second;
		rec_min = std::min(rec_value, rec_min);
	}

	if (rec_min < std::numeric_limits<ai_t>::max()) {
		return i_c + rec_min;
	} else {
		return std::numeric_limits<ai_t>::max();
	}
}
void update_startpoint_case_one_naive(case_one_index_naive &I, const anchor_t &a_j, const ai_t j_cost)
{
	// do nothing
}
void update_endpoint_case_one_naive(case_one_index_naive &I, const anchor_t &a_j, const ai_t j_cost)
{
	ai_t j_a = get<0>(a_j);
	ai_t j_c = get<1>(a_j);
	ai_t j_d = get<1>(a_j) + get<2>(a_j);
	ai_t j_diag = j_a - j_c;
	ai_t rec_value = j_cost - j_d;

	if (I.recursive_values.contains(j_diag)) {
		I.recursive_values[j_diag] = std::min(I.recursive_values[j_diag], rec_value);
	} else {
		I.recursive_values[j_diag] = rec_value;
	}
}

struct case_one_index init_case_one(ai_t Tlength, ai_t Qlength)
{
	return case_one_index({ MinSegmentTree<ai_t>(-Qlength, Tlength) });
}

ai_t compute_case_one(const case_one_index &I, const anchor_t &a_i)
{
	const ai_t i_a = get<0>(a_i);
	const ai_t i_c = get<1>(a_i);
	const ai_t i_diag = i_a - i_c;

	// naive range min query
	const ai_t rec_min = I.recursive_values.query(i_diag + 1, I.recursive_values.maxquery);

	if (rec_min < std::numeric_limits<ai_t>::max()) {
		return i_c + rec_min;
	} else {
		return std::numeric_limits<ai_t>::max();
	}
}

void update_startpoint_case_one(case_one_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	// do nothing
}

void update_endpoint_case_one(case_one_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	const ai_t j_a = get<0>(a_j);
	const ai_t j_c = get<1>(a_j);
	const ai_t j_d = get<1>(a_j) + get<2>(a_j);
	const ai_t j_diag = j_a - j_c;
	const ai_t rec_value = j_cost - j_d;

	I.recursive_values.update(j_diag, rec_value);
}

ai_t compute_case_one_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs)
{
	ai_t cost = std::numeric_limits<ai_t>::max();

	ai_t i_a = get<0>(anchors[i]);
	ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
	ai_t i_c = get<1>(anchors[i]);
	ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

	// anchor j < anchor i 
	for(ai_t j = i - 1; j > 0; j--) // dummy start anchor is handled separately
	{
		ai_t j_a = get<0>(anchors[j]);
		ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]);
		ai_t j_c = get<1>(anchors[j]);
		ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]);

		if (costs[i] < std::numeric_limits<ai_t>::max() and
		    j_b <= i_a and j_d <= i_c and j_a - j_c > i_a - i_c) // case 1
			cost = min(cost, costs[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d));
	}

	return cost;
}


} // namespace algo

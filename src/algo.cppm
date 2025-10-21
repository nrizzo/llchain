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

using std::vector;
using std::tuple, std::pair;
using std::max, std::min;
using std::cerr, std::endl;
using std::map;
using std::multiset;
using utils::anchor_t;

namespace algo {
/*
 * solves chainx-precedence colinear chaining via DP
 * NB: O(n^2) time
 * NB: assumes anchors contains dummies (see place_dummy_anchors)
 * NB: assumes sorted anchors (see sort_anchors)
 */
export void chainx_global_naive(
		const vector<anchor_t> &anchors,
		vector<long long> &costs_out,
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
		vector<long long> &costs_out,
		vector<anchor_t> &chain_out,
		const vector<long long> &correct_costs
	);

export void chainx_global_naive(const vector<anchor_t> &anchors, vector<long long> &costs_out, vector<anchor_t> &chain_out)
{
	long long n = anchors.size();
	costs_out = vector<long long>(n, 0);
	vector<long long> backtracks(n, -1);
	chain_out.clear();

	// compute optimal costs (ChainX precedence)
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

			if (costs_out[i] < std::numeric_limits<long long>::max() && i_a < j_a && i_b < j_b && i_c < j_c && i_d < j_d)
			{
				long long gap1 = max((long long)0, j_a - i_b - 1);
				long long gap2 = max((long long)0, j_c - i_d - 1);
				long long g = max(gap1,gap2);

				long long overlap1 = max((long long)0, i_b - j_a + 1);
				long long overlap2 = max((long long)0, i_d - j_c + 1);
				long long o = abs(overlap1 - overlap2);

				if (costs_out[i] + g + o < find_min_cost) backtrack = i;
				find_min_cost = min(find_min_cost, costs_out[i] + g + o);
			}
		}
		//save optimal cost at offset j
		costs_out[j] = find_min_cost;
		backtracks[j] = backtrack;
	}

	// return optimal chain
	for (long long i = n - 1; backtracks[i] >= 0; i = backtracks[i]) {
		chain_out.push_back(anchors[i]);
	}
	chain_out.push_back(anchors[0]);
	std::reverse(chain_out.begin(), chain_out.end());
}

struct case_one_index {
	// this 1D RmQ data structure needs to support updates to better values only
	// we know the position of all future queries beforehand, but not their values
	// TODO implement an efficient version
	std::map<long long, long long> recursive_values; // diagonal -> min(C[j] - q_e^j)
};

struct case_one_index init_case_one()
{
	return case_one_index();
}

long long compute_case_one(const case_one_index &I, const anchor_t &a_i)
{
	long long i_a = get<0>(a_i);
	long long i_c = get<1>(a_i);
	long long i_diag = i_a - i_c;

	// naive range min query
	long long rec_min = std::numeric_limits<long long>::max();
	for (auto it = I.recursive_values.upper_bound(i_diag); it != I.recursive_values.end(); it++) {
		const long long rec_value = it->second;
		rec_min = std::min(rec_value, rec_min);
	}

	if (rec_min < std::numeric_limits<long long>::max()) {
		return i_c + rec_min;
	} else {
		return std::numeric_limits<long long>::max();
	}
}

void update_startpoint_case_one(case_one_index &I, const anchor_t &a_j, const long long j_cost)
{
	// do nothing
}

void update_endpoint_case_one(case_one_index &I, const anchor_t &a_j, const long long j_cost)
{
	long long j_a = get<0>(a_j);
	long long j_c = get<1>(a_j);
	long long j_d = get<1>(a_j) + get<2>(a_j);
	long long j_diag = j_a - j_c;
	long long rec_value = j_cost - j_d;

	if (I.recursive_values.contains(j_diag)) {
		I.recursive_values[j_diag] = std::min(I.recursive_values[j_diag], rec_value);
	} else {
		I.recursive_values[j_diag] = rec_value;
	}
}

long long naive_compute_case_one(const vector<anchor_t> &anchors, long long i, const vector<long long> &costs_out)
{
	// naively compute what we THINK we are computing in case one 
	long long find_min_cost = std::numeric_limits<long long>::max();

	long long i_a = get<0>(anchors[i]);
	long long i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
	long long i_c = get<1>(anchors[i]);
	long long i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

	// anchor j < anchor i 
	for(long long j = i - 1; j > 0; j--) // dummy anchor[0] is treated separately!
	{
		long long j_a = get<0>(anchors[j]);
		long long j_b = get<0>(anchors[j]) + get<2>(anchors[j]) - 1;
		long long j_c = get<1>(anchors[j]);
		long long j_d = get<1>(anchors[j]) + get<2>(anchors[j]) - 1;

		if (costs_out[j] < std::numeric_limits<long long>::max() and 
				j_b < i_a and j_d < i_c and j_a - j_c > i_a - i_c) // case 1
		{
			long long gap1 = max((long long)0, i_a - j_b - 1);
			long long gap2 = max((long long)0, i_c - j_d - 1);
			long long g = max(gap1,gap2);

			long long overlap1 = max((long long)0, j_b - i_a + 1);
			long long overlap2 = max((long long)0, j_d - i_c + 1);
			long long o = abs(overlap1 - overlap2);

			find_min_cost = min(find_min_cost, costs_out[j] + g + o);
		}
	}

	return find_min_cost;
}

enum line
{
	horizontal,
	diagonal
};

struct case_two_index {
	// TODO implement an efficient version
	const vector<anchor_t> &anchors;
	vector<pair<long long,line>> delimiting_lines; // element is pair of anchor index and line type
	vector<long long> optimal_recursion_values;
};

struct case_two_index init_case_two(const vector<anchor_t> &anchors) {
	case_two_index I = {anchors, {{0, horizontal}, {anchors.size()-1,horizontal}}, {std::numeric_limits<long long>::max()}};
	return I;
}

long long get_c(case_two_index &I, long long line_index, long long current_a)
{
	assert(line_index >= 0 and line_index < I.delimiting_lines.size());
	long long i = get<0>(I.delimiting_lines[line_index]);
	if (get<1>(I.delimiting_lines[line_index]) == horizontal) {
		long long i_d = get<1>(I.anchors[i]) + get<2>(I.anchors[i]);
		return i_d;
	} else { // diagonal
		long long i_diag = get<0>(I.anchors[i]) - get<1>(I.anchors[i]);
		return current_a - i_diag;
	}
}

void prune_shadowed_delimiting_lines(case_two_index &I, long long sweeping_line_a)
{
	if (I.delimiting_lines.size() > 0) {
		assert(I.delimiting_lines.size() == I.optimal_recursion_values.size() + 1);
		assert(I.delimiting_lines[0].second == horizontal);
		for(long long l = 1; l < I.delimiting_lines.size() - 1; l++) { // l>=2?
			// has a line reached the following one?
			pair<long long,line> l1 = I.delimiting_lines[l];
			pair<long long,line> l2 = I.delimiting_lines[l + 1];
			if (l1.second == diagonal and l2.second == horizontal) {
				long long l1_diag = get<0>(I.anchors[l1.first]) - get<1>(I.anchors[l1.first]);
				long long l2_d = get<1>(I.anchors[l2.first]) + get<2>(I.anchors[l2.first]);

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

long long compute_case_two(case_two_index &I, const anchor_t &anchor_i)
{
	// TODO binary search
	long long i_a = get<0>(anchor_i);
	long long i_c = get<1>(anchor_i);
	prune_shadowed_delimiting_lines(I, i_a);

	// no recursion available
	if (I.delimiting_lines.size() == 0 or
			i_c < get_c(I, 0, i_a) or
			i_c > get_c(I, I.delimiting_lines.size() - 1, i_a))
		return std::numeric_limits<long long>::max();

	// recursion available
	long long i = 0;
	while (i_c >= get_c(I, i, i_a)) i++;
	if (I.delimiting_lines[i-1].second == diagonal and i_c == get_c(I, i-1, i_a)) i--;
	if (I.optimal_recursion_values[i-1] != std::numeric_limits<long long>::max()) {
		return I.optimal_recursion_values[i-1] + i_a;
	} else {
		return std::numeric_limits<long long>::max();
	}
}

void update_startpoint_case_two(case_two_index &I, const anchor_t &a_j, const long long j_cost)
{
	// do nothing
}

void update_endpoint_case_two(case_two_index &I, long long j, const anchor_t &a_j, const long long j_cost)
{
	long long j_b = get<0>(a_j) + get<2>(a_j);
	long long j_d = get<1>(a_j) + get<2>(a_j);
	prune_shadowed_delimiting_lines(I, j_b);

	// value of the two lines to potentially insert
	long long recursion_value = j_cost - j_b;

	// find position in array
	long long i = 0;
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

long long naive_compute_case_two(const vector<anchor_t> &anchors, long long i, const vector<long long> &costs_out)
{
	// naively compute what we THINK we are computing in case two
	long long find_min_cost = std::numeric_limits<long long>::max();

	long long i_a = get<0>(anchors[i]);
	long long i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
	long long i_c = get<1>(anchors[i]);
	long long i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

	// anchor j < anchor i 
	for(long long j = i - 1; j > 0; j--) // dummy anchor[0] is treated separately!
	{
		long long j_a = get<0>(anchors[j]);
		long long j_b = get<0>(anchors[j]) + get<2>(anchors[j]) - 1;
		long long j_c = get<1>(anchors[j]);
		long long j_d = get<1>(anchors[j]) + get<2>(anchors[j]) - 1;

		if (costs_out[j] < std::numeric_limits<long long>::max() and 
				j_b < i_a and j_d < i_c and j_a - j_c <= i_a - i_c) // case 2
		{
			long long gap1 = max((long long)0, i_a - j_b - 1);
			long long gap2 = max((long long)0, i_c - j_d - 1);
			long long g = max(gap1,gap2);

			long long overlap1 = max((long long)0, j_b - i_a + 1);
			long long overlap2 = max((long long)0, j_d - i_c + 1);
			long long o = abs(overlap1 - overlap2);

			find_min_cost = min(find_min_cost, costs_out[j] + g + o);
		}
	}

	return find_min_cost;
}

struct case_three_index {
	// TODO implement an efficient version
	std::map<long long, multiset<long long>> recursive_values; // diagonal -> multiset of values C[j] + diag_j
};

struct case_three_index init_case_three()
{
	return case_three_index();
}

long long compute_case_three(const case_three_index &I, const anchor_t &a_i)
{
	long long i_a = get<0>(a_i);
	long long i_c = get<1>(a_i);
	long long i_diag = i_a - i_c;

	// naive range min query
	long long rec_min = std::numeric_limits<long long>::max();
	for (auto it = I.recursive_values.upper_bound(i_diag); it != I.recursive_values.end(); it++) {
		const multiset<long long> rec_values = it->second;
		rec_min = min(*min_element(rec_values.begin(), rec_values.end()), rec_min);
	}

	if (rec_min < std::numeric_limits<long long>::max()) {
		return rec_min - i_diag;
	} else {
		return std::numeric_limits<long long>::max();
	}
}

void update_startpoint_case_three(case_three_index &I, const anchor_t &a_j, const long long j_cost)
{
	long long j_a = get<0>(a_j);
	long long j_c = get<1>(a_j);
	long long j_diag = j_a - j_c;
	long long rec_value = j_cost + j_diag;

	if (I.recursive_values.contains(j_diag)) {
		I.recursive_values[j_diag].insert(rec_value);
	} else {
		I.recursive_values[j_diag] = { rec_value };
	}
}

void update_endpoint_case_three(case_three_index &I, const anchor_t &a_j, const long long j_cost)
{
	long long j_a = get<0>(a_j);
	long long j_c = get<1>(a_j);
	long long j_diag = j_a - j_c;
	long long rec_value = j_cost + j_diag;

	I.recursive_values[j_diag].erase(
			I.recursive_values[j_diag].find(rec_value)
		);
	if (I.recursive_values[j_diag].size() == 0) {
		I.recursive_values.erase(j_diag);
	}
}

long long naive_compute_case_three(const vector<anchor_t> &anchors, long long i, const vector<long long> &costs_out)
{
	// naively compute what we THINK we are computing in case three
	long long find_min_cost = std::numeric_limits<long long>::max();

	long long i_a = get<0>(anchors[i]);
	long long i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
	long long i_c = get<1>(anchors[i]);
	long long i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

	// anchor j < anchor i 
	for(long long j = i - 1; j > 0; j--) // dummy anchor[0] is treated separately!
	{
		long long j_a = get<0>(anchors[j]);
		long long j_b = get<0>(anchors[j]) + get<2>(anchors[j]) - 1;
		long long j_c = get<1>(anchors[j]);
		long long j_d = get<1>(anchors[j]) + get<2>(anchors[j]) - 1;

		if (costs_out[i] < std::numeric_limits<long long>::max() and 
				j_a <= i_a and i_a <= j_b and j_a - j_c > i_a - i_c) // case 3
		{
			long long gap1 = max((long long)0, i_a - j_b - 1);
			long long gap2 = max((long long)0, i_c - j_d - 1);
			long long g = max(gap1,gap2);

			long long overlap1 = max((long long)0, j_b - i_a + 1);
			long long overlap2 = max((long long)0, j_d - i_c + 1);
			long long o = abs(overlap1 - overlap2);

			find_min_cost = min(find_min_cost, costs_out[j] + g + o);
		}
	}

	return find_min_cost;
}

struct case_four_index {
	// TODO implement an efficient version
	std::map<long long, multiset<long long>> recursive_values; // diagonal -> multiset of values C[j] + diag_j
};

struct case_four_index init_case_four()
{
	return case_four_index();
}

long long compute_case_four(const case_four_index &I, const anchor_t &a_i)
{
	long long i_a = get<0>(a_i);
	long long i_c = get<1>(a_i);
	long long i_diag = i_a - i_c;

	// naive range min query
	long long rec_min = std::numeric_limits<long long>::max();
	for (auto it = I.recursive_values.begin(); it != I.recursive_values.upper_bound(i_diag); it++) {
		const multiset<long long> rec_values = it->second;
		rec_min = min(*min_element(rec_values.begin(), rec_values.end()), rec_min);
	}

	if (rec_min < std::numeric_limits<long long>::max()) {
		return rec_min + i_diag;
	} else {
		return std::numeric_limits<long long>::max();
	}
}

void update_startpoint_case_four(case_four_index &I, const anchor_t &a_j, const long long j_cost)
{
	long long j_a = get<0>(a_j);
	long long j_c = get<1>(a_j);
	long long j_diag = j_a - j_c;
	long long rec_value = j_cost - j_diag;

	if (I.recursive_values.contains(j_diag)) {
		I.recursive_values[j_diag].insert(rec_value);
	} else {
		I.recursive_values[j_diag] = { rec_value };
	}
}

void update_endpoint_case_four(case_four_index &I, const anchor_t &a_j, const long long j_cost)
{
	long long j_a = get<0>(a_j);
	long long j_c = get<1>(a_j);
	long long j_diag = j_a - j_c;
	long long rec_value = j_cost - j_diag;

	I.recursive_values[j_diag].erase(
			I.recursive_values[j_diag].find(rec_value)
		);
	if (I.recursive_values[j_diag].size() == 0) {
		I.recursive_values.erase(j_diag);
	}
}

long long naive_compute_case_four(const vector<anchor_t> &anchors, long long i, const vector<long long> &costs_out)
{
	// naively compute what we THINK we are computing in case four
	long long find_min_cost = std::numeric_limits<long long>::max();

	long long i_a = get<0>(anchors[i]);
	long long i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
	long long i_c = get<1>(anchors[i]);
	long long i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

	// anchor j < anchor i 
	for(long long j = i - 1; j > 0; j--) // dummy anchor[0] is treated separately!
	{
		long long j_a = get<0>(anchors[j]);
		long long j_b = get<0>(anchors[j]) + get<2>(anchors[j]) - 1;
		long long j_c = get<1>(anchors[j]);
		long long j_d = get<1>(anchors[j]) + get<2>(anchors[j]) - 1;

		if (costs_out[i] < std::numeric_limits<long long>::max() and 
				j_a <= i_a and i_a <= j_b and j_a - j_c <= i_a - i_c) // case 4
		{
			long long gap1 = max((long long)0, i_a - j_b - 1);
			long long gap2 = max((long long)0, i_c - j_d - 1);
			long long g = max(gap1,gap2);

			long long overlap1 = max((long long)0, j_b - i_a + 1);
			long long overlap2 = max((long long)0, j_d - i_c + 1);
			long long o = abs(overlap1 - overlap2);

			find_min_cost = min(find_min_cost, costs_out[j] + g + o);
		}
	}

	return find_min_cost;
}

struct case_five_index {
	// TODO implement an efficient version
	std::map<long long, long long> recursive_values; // j_d -> min(C[j] - diag_j)
};

struct case_five_index init_case_five()
{
	return case_five_index();
}

long long compute_case_five(const case_five_index &I, const anchor_t &a_i) {
	long long i_a = get<0>(a_i);
	long long i_c = get<1>(a_i);
	long long i_d = get<1>(a_i) + get<2>(a_i);
	long long i_diag = i_a - i_c;

	// naive range min query
	long long rec_min = std::numeric_limits<long long>::max();
	for (auto it = I.recursive_values.lower_bound(i_c); it != I.recursive_values.upper_bound(i_d - 1); it++) {
		const long long rec_value = it->second;
		rec_min = std::min(rec_value, rec_min);
	}

	if (rec_min < std::numeric_limits<long long>::max()) {
		return rec_min + i_diag;
	} else {
		return std::numeric_limits<long long>::max();
	}
}

void update_startpoint_case_five(case_five_index &I, const anchor_t &a_j, const long long j_cost)
{
	// do nothing
}

void update_endpoint_case_five(case_five_index &I, const anchor_t &a_j, const long long j_cost)
{
	long long j_a = get<0>(a_j);
	long long j_c = get<1>(a_j);
	long long j_d = get<1>(a_j) + get<2>(a_j);
	long long j_diag = j_a - j_c;
	long long rec_value = j_cost - j_diag;

	if (I.recursive_values.contains(j_d-1)) {
		I.recursive_values[j_d-1] = std::min(I.recursive_values[j_d-1], rec_value);
	} else {
		I.recursive_values[j_d-1] = rec_value;
	}
}

long long naive_compute_case_five(const vector<anchor_t> &anchors, long long i, const vector<long long> &costs_out)
{
	//compute cost[i] here
	long long find_min_cost = std::numeric_limits<long long>::max();

	long long i_a = get<0>(anchors[i]);
	long long i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
	long long i_c = get<1>(anchors[i]);
	long long i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

	// anchor j < anchor i 
	for(long long j = i - 1; j > 0; j--) // dummy anchor[0] is treated separately!
	{
		long long j_a = get<0>(anchors[j]);
		long long j_b = get<0>(anchors[j]) + get<2>(anchors[j]) - 1;
		long long j_c = get<1>(anchors[j]);
		long long j_d = get<1>(anchors[j]) + get<2>(anchors[j]) - 1;

		if (costs_out[i] < std::numeric_limits<long long>::max() and 
				//j_a < i_a and j_b < i_b and j_c < i_c and j_d < i_d and // TODO remove this?
				j_b < i_a and i_c <= j_d and j_d <= i_d) // case 5
		{
			long long gap1 = max((long long)0, i_a - j_b - 1);
			long long gap2 = max((long long)0, i_c - j_d - 1);
			long long g = max(gap1,gap2);

			long long overlap1 = max((long long)0, j_b - i_a + 1);
			long long overlap2 = max((long long)0, j_d - i_c + 1);
			long long o = abs(overlap1 - overlap2);

			find_min_cost = min(find_min_cost, costs_out[j] + g + o);
		}
	}

	return find_min_cost;
}

export void solve_global_linearithmic_naive(
		const vector<anchor_t> &anchors,
		vector<long long> &costs_out,
		vector<anchor_t> &chain_out,
		const vector<long long> &correct_costs
	)
{
	// NB anchors down here are closed-open intervals [i_a..i_b), [i_c..i_d) MOST OF THE TIME TODO
	long long n = anchors.size();
	costs_out = vector<long long>(n, 0);
	chain_out.clear();

	vector<long long> points; // horizontal line sweep (+j means start, -j means end)
	points.reserve(2 * n - 2);
	for (long long i = 1; i < n; i++) // skip starting dummy anchor
		points.push_back(-i);
	for (long long i = 1; i < n; i++) // skip starting dummy anchor
		points.push_back(i);
	std::stable_sort(points.begin(), points.end(),
			[&](const long long i,
				const long long j) -> bool
			{
			return (((i >= 0) ? std::get<0>(anchors[i]) : std::get<0>(anchors[-i]) + std::get<2>(anchors[-i])) <
					((j >= 0) ? std::get<0>(anchors[j]) : std::get<0>(anchors[-j]) + std::get<2>(anchors[-j])));
			});

	case_one_index   I_one   = init_case_one();
	case_two_index   I_two   = init_case_two(anchors);
	case_three_index I_three = init_case_three();
	case_four_index  I_four  = init_case_four();
	case_five_index  I_five  = init_case_five();

	// compute optimal costs (ChainX RELAXED precedence)
	for (long long point : points) {
		if (point >= 0) { // startpoint
			long long i = point;
			long long i_a = get<0>(anchors[i]);
			long long i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
			long long i_c = get<1>(anchors[i]);
			long long i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

			// compute cost[i]
			long long cost = costs_out[0] + max(i_a, i_c); // connect(a_0, a_i)
			assert(compute_case_one  (I_one,   anchors[i]) == naive_compute_case_one  (anchors, i, costs_out));
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
			long long i = -point;
			update_endpoint_case_one  (I_one,    anchors[i], costs_out[i]);
			update_endpoint_case_two  (I_two, i, anchors[i], costs_out[i]);
			update_endpoint_case_three(I_three,  anchors[i], costs_out[i]);
			update_endpoint_case_four (I_four,   anchors[i], costs_out[i]);
			update_endpoint_case_five (I_five,   anchors[i], costs_out[i]);
		}
	}
}
} // namespace algo

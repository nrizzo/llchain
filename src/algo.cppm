export module algo;

import <vector>;
import <tuple>;
import <algorithm>; // std::sort
import <iostream>;
import <map>;
import <set>;
import <list>;
import <cassert>;
import utils;
import MinSegmentTree;

using std::vector;
using std::tuple, std::pair;
using std::max, std::min, std::abs;
using std::cerr, std::endl;
using std::map, std::multimap;
using std::set, std::multiset;
using std::list;
using utils::anchor_t, utils::connect, utils::connect_Qgap, utils::chainx_precedes, utils::weak_precedes;
typedef utils::anchor_index_t ai_t;

namespace algo {
export enum chaining_mode { global, semiglobal };
/*
 * solves chainx-precedence colinear chaining via DP, returns an optimal chain
 *   via backtracking; adapted from github.com/at-cg/ChainX
 * NB: O(n^2) time
 * NB: assumes anchors contains dummies (see place_dummy_anchors)
 * NB: assumes sorted anchors (see sort_anchors)
 */
export void chainx_naive(
		const vector<anchor_t> &anchors,
		const chaining_mode m,
		vector<ai_t> &costs_out,
		vector<anchor_t> &chain_out
	);

/*
 * solves (relaxed chainx-precedence) colinear chaining via DP and mimicks the
 *   solution/data structures of our linearithmic-time prototype
 * NB: assumes anchors contains dummies (see place_dummy_anchors)
 * NB: assumes sorted anchors (see sort_anchors)
 * NB: assumes there are no perfect chains between the anchors
 */
export void solve_linearithmic(
		const vector<anchor_t> &anchors,
		const ai_t Tlength,
		const ai_t Qlength,
		const chaining_mode m,
		vector<ai_t> &costs_out
	);
/*
 * debug version that checks that all partial costs are correct
 * NB: O(n^2) time or worse
 */
export void solve_linearithmic_debug(
		const vector<anchor_t> &anchors,
		const ai_t Tlength,
		const ai_t Qlength,
		const chaining_mode m,
		vector<ai_t> &costs_out,
		const vector<ai_t> &correct_costs
	);

/*
 * backtrack an optimal chainx/weak chain (including the dummy anchors), given
 *   the partial chaining costs
 * NB: assumes anchors contains dummies (see place_dummy_anchors)
 * NB: assumes sorted anchors (see sort_anchors)
 * NB: assumes the DP costs have been found with the corresponding precedence
 */
export void chainx_backtrack(
		const vector<anchor_t> &anchors,
		const vector<ai_t> &costs,
		const chaining_mode m,
		vector<anchor_t> &chain_out
	);
export void weak_backtrack(
		const vector<anchor_t> &anchors,
		const vector<ai_t> &costs,
		const chaining_mode m,
		vector<anchor_t> &chain_out
	);

/*
 * computes cost of a chain
 * NB: assumes anchors contains dummies (see place_dummy_anchors)
 * NB: assumes the chain is colinear (see chainx_precedes and weak_precedes)
 */
export ai_t compute_chain_cost(const vector<anchor_t> &chain, const chaining_mode m);

// start of implementation
export
void chainx_naive(
		const vector<anchor_t> &anchors,
		const chaining_mode m,
		vector<ai_t> &costs_out,
		vector<anchor_t> &chain_out
) {
	const ai_t n = anchors.size();
	costs_out = vector<ai_t>(n); // partial chainx-prec DP costs
	vector<ai_t> backtracks(n, -1); // index of an optimal previous anchor
	chain_out.clear();

	costs_out[0] = 0; // first dummy anchor
	for (ai_t i = 1; i < ((m == global) ? n : n-1); i++) {
		ai_t i_cost = std::numeric_limits<ai_t>::max();
		ai_t backtrack = -1;

		const ai_t i_a = get<0>(anchors[i]);
		const ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
		const ai_t i_c = get<1>(anchors[i]);
		const ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

		if (m == global) {
			i_cost = connect(anchors[0], i_a, i_b, i_c, i_d);
		} else {
			i_cost = connect_Qgap(anchors[0], i_c);
		}
		backtrack = 0;

		for (ai_t j = i - 1; j > 0; j--) { // anchor j < anchor i, first dummy anchor is treated separately
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
	if (m == semiglobal) {
		ai_t final_cost = std::numeric_limits<ai_t>::max();
		ai_t backtrack = -1;
		const ai_t final_c = get<1>(anchors[n-1]);
		for (ai_t j = 1; j < n - 1; j++) {
			if (costs_out[j] < std::numeric_limits<ai_t>::max()) {
				const ai_t c = costs_out[j] + connect_Qgap(anchors[j], final_c);
				if (c < final_cost) backtrack = j;
				final_cost = min(final_cost, c);
			}
		}
		costs_out[n-1] = final_cost;
		backtracks[n-1] = backtrack;
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
 * helper function that considers all recursive cases considered by case one
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

/*
 * case two: gap-gap case, bigger or equal-length gap in T
 *   requires a complex data structure that handles all boundaries obtained
 *   by propagating case 2 recursions to the right (horizontal or diagonal
 *   lines)
 *   TODO implement backtracking
 */
enum line_t { horizontal, diagonal };
typedef tuple<ai_t, line_t, ai_t> l_t; // (anchor index j, hor/diag, C[k] - t_e^k)
struct case_two_index {
	const vector<anchor_t>& anchors;
	list<l_t> delimiting_lines; // top-to-bottom ordered list of lines/rec values involved in case 2
	vector<list<l_t>::iterator> anchor_diag_to_list; // [j] contains pointer to diagonal   list element of anchor[j]
	vector<list<l_t>::iterator> anchor_hor_to_list;  // [j] contains pointer to horizontal list element of anchor[j]
	multimap<ai_t,ai_t> updates; // event T-position -> anchor j whose diagonal line might cross the next line
	map<ai_t, ai_t> active_horizontal_lines; // horizontal line -> anchor j
	map<ai_t, ai_t> active_diagonal_lines;   // diagonal line -> anchor j
};
struct case_two_index init_case_two(const vector<anchor_t> &anchors);
ai_t compute_case_two(case_two_index &I, const anchor_t &anchor_i);
void update_startpoint_case_two(case_two_index &I, const anchor_t &a_j, const ai_t j_cost);
void  update_endpoint_case_two (case_two_index &I, ai_t j, const anchor_t &a_j, const ai_t j_cost);

struct case_two_index_naive {
	const vector<anchor_t>& anchors;
	vector<pair<ai_t,line_t>> delimiting_lines;
	vector<ai_t> optimal_recursion_values;
};
struct case_two_index_naive init_case_two_naive(const vector<anchor_t> &anchors);
ai_t compute_case_two_naive(case_two_index_naive &I, const anchor_t &anchor_i);
void update_startpoint_case_two_naive(case_two_index_naive &I, const anchor_t &a_j, const ai_t j_cost);
void  update_endpoint_case_two_naive (case_two_index_naive &I, ai_t j, const anchor_t &a_j, const ai_t j_cost);

/*
 * helper function that considers all recursive cases considered by case two
 */
ai_t compute_case_two_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs_out);

/*
 * helper function that gets the Q-position of a horizontal or diagonal line contained in the index
 *   for case 2
 */
ai_t get_c(case_two_index_naive &I, ai_t line_index, ai_t current_a);

/*
 * helper function to handle intersection events for case 2 and to remove lines
 *   that have been crossed
 * NB: the naive version can take O(n^2) per call (but is probably O(n)
 *     in practice)
 */
void prune_shadowed_delimiting_lines_naive(case_two_index_naive &I, ai_t sweeping_line_a);
void prune_shadowed_delimiting_lines(case_two_index &I, const ai_t sweeping_line_a);

/*
 * case three: overlap in T, bigger diagonal
 *   requires a 1-dimensional range minimum query data structure where the
 *   (static) positions are anchor diagonals and the values (dynamic, better
 *   values or removal ops only) are positive
 */
struct case_three_index {
	MinSegmentTree<ai_t> recursive_values; // diagonal -> current C[j] + diag_j
};
struct case_three_index init_case_three(ai_t Tlength, ai_t Qlength);
ai_t compute_case_three(const case_three_index &I, const anchor_t &a_i);
void update_startpoint_case_three(case_three_index &I, const anchor_t &a_j, const ai_t j_cost);
void  update_endpoint_case_three (case_three_index &I, const anchor_t &a_j, const ai_t j_cost);

/*
 * helper function that considers all recursive cases considered by case three
 */
ai_t compute_case_three_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs);

// unused functions below
struct case_three_index_naive {
	std::map<ai_t, multiset<ai_t>> recursive_values; // diagonal -> values C[j] + diag_j
};
struct case_three_index_naive init_case_three_naive();
ai_t compute_case_three_naive(const case_three_index_naive &I, const anchor_t &a_i);
void update_startpoint_case_three_naive(case_three_index_naive &I, const anchor_t &a_j, const ai_t j_cost);
void update_endpoint_case_three_naive(case_three_index_naive &I, const anchor_t &a_j, const ai_t j_cost);

/*
 * case four: overlap in T, smaller diagonal
 *   requires a 1-dimensional range minimum query data structure where the
 *   (static) positions are anchor diagonals and the values (dynamic, better
 *   values or removal ops only) can be negative
 */
struct case_four_index {
	MinSegmentTree<ai_t> recursive_values; // diagonal -> current C[j] - diag_j
};
struct case_four_index init_case_four(ai_t Tlength, ai_t Qlength);
ai_t compute_case_four(const case_four_index &I, const anchor_t &a_i);
void update_startpoint_case_four(case_four_index &I, const anchor_t &a_j, const ai_t j_cost);
void  update_endpoint_case_four (case_four_index &I, const anchor_t &a_j, const ai_t j_cost);

/*
 * helper function that considers all recursive cases considered by case four
 */
ai_t compute_case_four_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs_out);

// unused functions below
struct case_four_index_naive {
	std::map<ai_t, multiset<ai_t>> recursive_values; // diagonal -> multiset of values C[j] + diag_j
};
struct case_four_index_naive init_case_four_naive();
ai_t compute_case_four_naive(const case_four_index_naive &I, const anchor_t &a_i);
void update_startpoint_case_four_naive(case_four_index_naive &I, const anchor_t &a_j, const ai_t j_cost);
void update_endpoint_case_four_naive(case_four_index_naive &I, const anchor_t &a_j, const ai_t j_cost);

/*
 * case five: gap in T, overlap in Q
 *   requires a 1-dimensional range minimum query data structure where the
 *   (static) positions are Q positions and the values (semi-dynamic, better
 *   values only) can be negative
 */
struct case_five_index {
	MinSegmentTree<ai_t> recursive_values; // Q endpoint -> C[j] - diag_j
};
struct case_five_index init_case_five(ai_t Tlength, ai_t Qlength);
ai_t compute_case_five(const case_five_index &I, const anchor_t &a_i);
void update_startpoint_case_five(case_five_index &I, const anchor_t &a_j, const ai_t j_cost);
void  update_endpoint_case_five (case_five_index &I, const anchor_t &a_j, const ai_t j_cost);

/*
 * helper function that considers all recursive cases considered by case five
 */
ai_t compute_case_five_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs_out);

// unused functions below
struct case_five_index_naive {
	std::map<ai_t, ai_t> recursive_values; // j_d -> min(C[j] - diag_j)
};
struct case_five_index_naive init_case_five_naive();
ai_t compute_case_five_naive(const case_five_index_naive &I, const anchor_t &a_i);
void update_startpoint_case_five_naive(case_five_index_naive &I, const anchor_t &a_j, const ai_t j_cost);
void update_endpoint_case_five_naive(case_five_index_naive &I, const anchor_t &a_j, const ai_t j_cost);

export
void solve_linearithmic(
		const vector<anchor_t> &anchors,
		const ai_t Tlength,
		const ai_t Qlength,
		const chaining_mode m,
		vector<ai_t> &costs_out,
		vector<anchor_t> &chain_out
) {
	const ai_t n = anchors.size();
	costs_out = vector<ai_t>(n, 0);
	chain_out.clear();

	vector<ai_t> points; // horizontal line sweep (+i means T-startpoint of i-th anchor, -i means T-endpoint)
	points.reserve(2 * n - 3);
	for (ai_t i = 1; i < n - 1; i++) // skip both dummy anchors
		points.push_back(-i);
	for (ai_t i = 1; i < ((m == global) ? n : n-1); i++) // skip starting dummy anchor (and end anchor if in semi-global mode)
		points.push_back(i);
	std::stable_sort(points.begin(), points.end(),
			[&anchors](const ai_t i,
				const ai_t j) -> bool
			{
			return (((i >= 0) ? std::get<0>(anchors[i]) : std::get<0>(anchors[-i]) + std::get<2>(anchors[-i])) <
					((j >= 0) ? std::get<0>(anchors[j]) : std::get<0>(anchors[-j]) + std::get<2>(anchors[-j])));
			});

	case_one_index   I_one   = init_case_one(Tlength, Qlength);
	case_two_index   I_two   = init_case_two(anchors);
	case_three_index I_three = init_case_three(Tlength, Qlength);
	case_four_index  I_four  = init_case_four(Tlength, Qlength);
	case_five_index  I_five  = init_case_five(Tlength, Qlength);

	for (ai_t point : points) {
		if (point >= 0) { // startpoint
			ai_t i = point;
			ai_t i_a = get<0>(anchors[i]);
			ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
			ai_t i_c = get<1>(anchors[i]);
			ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

			// compute cost[i]
			ai_t cost = std::numeric_limits<ai_t>::max();
			if (m == global) {
				cost = connect(anchors[0], i_a, i_b, i_c, i_d);
			} else {
				cost = connect_Qgap(anchors[0], i_c);
			}

			cost = std::min(cost, compute_case_one  (I_one,   anchors[i]));
			cost = std::min(cost, compute_case_two  (I_two,   anchors[i]));
			cost = std::min(cost, compute_case_three(I_three, anchors[i]));
			cost = std::min(cost, compute_case_four (I_four,  anchors[i]));
			cost = std::min(cost, compute_case_five (I_five,  anchors[i]));
			costs_out[i] = cost;

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
	if (m == semiglobal) {
		ai_t final_cost = std::numeric_limits<ai_t>::max();
		ai_t backtrack = -1;
		const ai_t final_c = get<1>(anchors[n-1]);
		for (ai_t j = 1; j < n - 1; j++) {
			if (costs_out[j] < std::numeric_limits<ai_t>::max()) {
				const ai_t c = costs_out[j] + connect_Qgap(anchors[j], final_c);
				if (c < final_cost) backtrack = j;
				final_cost = min(final_cost, c);
			}
		}
		costs_out[n-1] = final_cost;
	}
}

export
void solve_linearithmic_debug(
		const vector<anchor_t> &anchors,
		const ai_t Tlength,
		const ai_t Qlength,
		const chaining_mode m,
		vector<ai_t> &costs_out,
		const vector<ai_t> &correct_costs
) {
	ai_t n = anchors.size();
	costs_out = vector<ai_t>(n, 0);

	vector<ai_t> points; // horizontal line sweep (+i means T-startpoint of i-th anchor, -i means T-endpoint)
	points.reserve(2 * n - 2);
	for (ai_t i = 1; i < n - 1; i++) // skip both dummy anchors
		points.push_back(-i);
	for (ai_t i = 1; i < ((m == global) ? n : n-1); i++) // skip starting dummy anchor (and end anchor if in semi-global mode)
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
	case_two_index_naive   I_two_naive   = init_case_two_naive(anchors);
	case_three_index I_three = init_case_three(Tlength, Qlength);
	case_four_index  I_four  = init_case_four(Tlength, Qlength);
	case_five_index  I_five  = init_case_five(Tlength, Qlength);

	for (ai_t point : points) {
		if (point >= 0) { // startpoint
			ai_t i = point;
			ai_t i_a = get<0>(anchors[i]);
			ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
			ai_t i_c = get<1>(anchors[i]);
			ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

			// compute cost[i]
			ai_t cost = std::numeric_limits<ai_t>::max();
			if (m == global) {
				cost = connect(anchors[0], i_a, i_b, i_c, i_d);
			} else {
				cost = connect_Qgap(anchors[0], i_c);
			}
			assert(compute_case_one  (I_one,   anchors[i]) == compute_case_one_debug (anchors, i, costs_out));
			assert(compute_case_two  (I_two,   anchors[i]) == compute_case_two_debug  (anchors, i, costs_out));
			assert(compute_case_two_naive  (I_two_naive,   anchors[i]) == compute_case_two_debug  (anchors, i, costs_out));
			assert(compute_case_three(I_three, anchors[i]) == compute_case_three_debug(anchors, i, costs_out));
			assert(compute_case_four (I_four,  anchors[i]) == compute_case_four_debug (anchors, i, costs_out));
			assert(compute_case_five (I_five,  anchors[i]) == compute_case_five_debug (anchors, i, costs_out));

			cost = std::min(cost, compute_case_one  (I_one,   anchors[i]));
			cost = std::min(cost, compute_case_two  (I_two,   anchors[i]));
			cost = std::min(cost, compute_case_three(I_three, anchors[i]));
			cost = std::min(cost, compute_case_four (I_four,  anchors[i]));
			cost = std::min(cost, compute_case_five (I_five,  anchors[i]));
			costs_out[i] = cost;
			assert(costs_out[i] <= correct_costs[i]);

			update_startpoint_case_one  (I_one,   anchors[i], costs_out[i]);
			update_startpoint_case_two  (I_two,   anchors[i], costs_out[i]);
			update_startpoint_case_two_naive  (I_two_naive,   anchors[i], costs_out[i]);
			update_startpoint_case_three(I_three, anchors[i], costs_out[i]);
			update_startpoint_case_four (I_four,  anchors[i], costs_out[i]);
			update_startpoint_case_five (I_five,  anchors[i], costs_out[i]);
		} else { // endpoint
			ai_t i = -point;
			update_endpoint_case_one  (I_one,    anchors[i], costs_out[i]);
			update_endpoint_case_two  (I_two, i, anchors[i], costs_out[i]);
			update_endpoint_case_two_naive  (I_two_naive, i, anchors[i], costs_out[i]);
			update_endpoint_case_three(I_three,  anchors[i], costs_out[i]);
			update_endpoint_case_four (I_four,   anchors[i], costs_out[i]);
			update_endpoint_case_five (I_five,   anchors[i], costs_out[i]);
		}
	}
	if (m == semiglobal) {
		ai_t final_cost = std::numeric_limits<ai_t>::max();
		const ai_t final_c = get<1>(anchors[n-1]);
		for (ai_t j = 1; j < n - 1; j++) {
			if (costs_out[j] < std::numeric_limits<ai_t>::max()) {
				const ai_t c = costs_out[j] + connect_Qgap(anchors[j], final_c);
				final_cost = min(final_cost, c);
			}
		}
		costs_out[n-1] = final_cost;
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

	const ai_t i_a = get<0>(anchors[i]);
	const ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
	const ai_t i_c = get<1>(anchors[i]);
	const ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

	// anchor j < anchor i
	for(ai_t j = i - 1; j > 0; j--) { // dummy start anchor is handled separately
		const ai_t j_a = get<0>(anchors[j]);
		const ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]);
		const ai_t j_c = get<1>(anchors[j]);
		const ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]);

		if (costs[j] < std::numeric_limits<ai_t>::max() and
		    j_b <= i_a and j_d <= i_c and j_a - j_c > i_a - i_c) // case 1
			cost = min(cost, costs[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d));
	}

	return cost;
}

struct case_two_index_naive init_case_two_naive(const vector<anchor_t> &anchors) {
	case_two_index_naive I = {anchors, {{0, horizontal}, {anchors.size()-1,horizontal}}, {std::numeric_limits<ai_t>::max()}};
	return I;
}

ai_t get_c(case_two_index_naive &I, ai_t line_index, ai_t current_a)
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

void prune_shadowed_delimiting_lines_naive(case_two_index_naive &I, ai_t sweeping_line_a)
{
	if (I.delimiting_lines.size() > 0) {
		assert(I.delimiting_lines.size() == I.optimal_recursion_values.size() + 1);
		assert(I.delimiting_lines[0].second == horizontal);
		for(ai_t l = 1; l < I.delimiting_lines.size() - 1; l++) { // l>=2?
			// has a line reached the following one?
			pair<ai_t,line_t> l1 = I.delimiting_lines[l];
			pair<ai_t,line_t> l2 = I.delimiting_lines[l + 1];
			if (l1.second == diagonal and l2.second == horizontal) {
				ai_t l1_diag = get<0>(I.anchors[l1.first]) - get<1>(I.anchors[l1.first]);
				ai_t l2_d = get<1>(I.anchors[l2.first]) + get<2>(I.anchors[l2.first]);

				if (sweeping_line_a - l1_diag >= l2_d) {
					if (l == I.delimiting_lines.size() - 2 or I.optimal_recursion_values[l-1] <= I.optimal_recursion_values[l+1]) {
						// diagonal shadows horizontal
						I.delimiting_lines.erase(I.delimiting_lines.begin() + l + 1);
						I.optimal_recursion_values.erase(I.optimal_recursion_values.begin() + l);
						l -= 1;
					} else {
						// horizontal shadows diagonal
						I.delimiting_lines.erase(I.delimiting_lines.begin() + l);
						I.optimal_recursion_values.erase(I.optimal_recursion_values.begin() + l);
						l = max(l - 2, (ai_t)0);
					}
				}
			}
		}
	}
}

ai_t compute_case_two_naive(case_two_index_naive &I, const anchor_t &anchor_i)
{
	ai_t i_a = get<0>(anchor_i);
	ai_t i_c = get<1>(anchor_i);
	prune_shadowed_delimiting_lines_naive(I, i_a);

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

void update_startpoint_case_two_naive(case_two_index_naive &I, const anchor_t &a_j, const ai_t j_cost)
{
	// do nothing
}

void update_endpoint_case_two_naive(case_two_index_naive &I, ai_t j, const anchor_t &a_j, const ai_t j_cost)
{
	ai_t j_b = get<0>(a_j) + get<2>(a_j);
	ai_t j_d = get<1>(a_j) + get<2>(a_j);
	prune_shadowed_delimiting_lines_naive(I, j_b);

	// value of the two lines to potentially insert
	ai_t recursion_value = j_cost - j_b;

	// find position in array
	ai_t i = 0;
	while (i < I.delimiting_lines.size() and j_d > get_c(I, i, j_b)) i++;

	if (I.delimiting_lines.size() == 0) {
		// TODO remove this branch, as it never gets called
		cerr << "DEBUG: does case 0 happen?" << endl;
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

ai_t compute_case_two_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs)
{
	ai_t cost = std::numeric_limits<ai_t>::max();

	ai_t i_a = get<0>(anchors[i]);
	ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
	ai_t i_c = get<1>(anchors[i]);
	ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

	// anchor j < anchor i
	for(ai_t j = i - 1; j > 0; j--) { // dummy start anchor is handled separately
		ai_t j_a = get<0>(anchors[j]);
		ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]);
		ai_t j_c = get<1>(anchors[j]);
		ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]);

		if (costs[j] < std::numeric_limits<ai_t>::max() and
		    j_b <= i_a and j_d <= i_c and j_a - j_c <= i_a - i_c) // case 2
			cost = min(cost, costs[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d));
	}

	return cost;
}

struct case_three_index_naive init_case_three_naive()
{
	return case_three_index_naive();
}

ai_t compute_case_three_naive(const case_three_index_naive &I, const anchor_t &a_i)
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

void update_startpoint_case_three_naive(case_three_index_naive &I, const anchor_t &a_j, const ai_t j_cost)
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

void update_endpoint_case_three_naive(case_three_index_naive &I, const anchor_t &a_j, const ai_t j_cost)
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

ai_t compute_case_three_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs)
{
	ai_t cost = std::numeric_limits<ai_t>::max();

	const ai_t i_a = get<0>(anchors[i]);
	const ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
	const ai_t i_c = get<1>(anchors[i]);
	const ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

	// anchor j < anchor i
	for(ai_t j = i - 1; j > 0; j--) { // dummy start anchor is handled separately
		const ai_t j_a = get<0>(anchors[j]);
		const ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]);
		const ai_t j_c = get<1>(anchors[j]);
		const ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]);

		if (costs[j] < std::numeric_limits<ai_t>::max() and
		    j_a <= i_a and i_a < j_b and j_a - j_c > i_a - i_c) // case 3
			cost = min(cost, costs[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d));
	}

	return cost;
}

struct case_three_index init_case_three(ai_t Tlength, ai_t Qlength)
{
	return case_three_index({ MinSegmentTree<ai_t>(-Qlength, Tlength) });
}

ai_t compute_case_three(const case_three_index &I, const anchor_t &a_i)
{
	const ai_t i_a = get<0>(a_i);
	const ai_t i_c = get<1>(a_i);
	const ai_t i_diag = i_a - i_c;

	const ai_t rec_min = I.recursive_values.query(i_diag + 1, I.recursive_values.maxquery);

	if (rec_min < std::numeric_limits<ai_t>::max()) {
		return rec_min - i_diag;
	} else {
		return std::numeric_limits<ai_t>::max();
	}
}

void update_startpoint_case_three(case_three_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	const ai_t j_a = get<0>(a_j);
	const ai_t j_c = get<1>(a_j);
	const ai_t j_diag = j_a - j_c;

	assert(I.recursive_values.query(j_diag, j_diag) == std::numeric_limits<ai_t>::max());
	I.recursive_values.update(j_diag, j_cost + j_diag);
}

void  update_endpoint_case_three (case_three_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	const ai_t j_a = get<0>(a_j);
	const ai_t j_c = get<1>(a_j);
	const ai_t j_diag = j_a - j_c;

	assert(I.recursive_values.query(j_diag, j_diag) == j_cost + j_diag);
	I.recursive_values.remove(j_diag);
}

struct case_four_index_naive init_case_four_naive()
{
	return case_four_index_naive();
}

ai_t compute_case_four_naive(const case_four_index_naive &I, const anchor_t &a_i)
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

void update_startpoint_case_four_naive(case_four_index_naive &I, const anchor_t &a_j, const ai_t j_cost)
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

void update_endpoint_case_four_naive(case_four_index_naive &I, const anchor_t &a_j, const ai_t j_cost)
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

ai_t compute_case_four_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs)
{
	ai_t cost = std::numeric_limits<ai_t>::max();

	const ai_t i_a = get<0>(anchors[i]);
	const ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
	const ai_t i_c = get<1>(anchors[i]);
	const ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

	// anchor j < anchor i
	for(ai_t j = i - 1; j > 0; j--) { // dummy start anchor is handled separately
		const ai_t j_a = get<0>(anchors[j]);
		const ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]);
		const ai_t j_c = get<1>(anchors[j]);
		const ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]);

		if (costs[j] < std::numeric_limits<ai_t>::max() and
		    j_a <= i_a and i_a < j_b and j_a - j_c < i_a - i_c) // case 4
			cost = min(cost, costs[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d));
	}

	return cost;
}

struct case_four_index init_case_four(ai_t Tlength, ai_t Qlength)
{
	return case_four_index({ MinSegmentTree<ai_t>(-Qlength, Tlength) });
}

ai_t compute_case_four(const case_four_index &I, const anchor_t &a_i)
{
	const ai_t i_a = get<0>(a_i);
	const ai_t i_c = get<1>(a_i);
	const ai_t i_diag = i_a - i_c;

	const ai_t rec_min = I.recursive_values.query(I.recursive_values.minquery, i_diag - 1);

	if (rec_min < std::numeric_limits<ai_t>::max()) {
		return rec_min + i_diag;
	} else {
		return std::numeric_limits<ai_t>::max();
	}
}

void update_startpoint_case_four(case_four_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	const ai_t j_a = get<0>(a_j);
	const ai_t j_c = get<1>(a_j);
	const ai_t j_diag = j_a - j_c;

	assert(I.recursive_values.query(j_diag, j_diag) == std::numeric_limits<ai_t>::max());
	I.recursive_values.update(j_diag, j_cost - j_diag);
}

void  update_endpoint_case_four (case_four_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	const ai_t j_a = get<0>(a_j);
	const ai_t j_c = get<1>(a_j);
	const ai_t j_diag = j_a - j_c;

	assert(I.recursive_values.query(j_diag, j_diag) == j_cost - j_diag);
	I.recursive_values.remove(j_diag);
}

struct case_five_index_naive init_case_five_naive()
{
	return case_five_index_naive();
}

ai_t compute_case_five_naive(const case_five_index_naive &I, const anchor_t &a_i)
{
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

void update_startpoint_case_five_naive(case_five_index_naive &I, const anchor_t &a_j, const ai_t j_cost)
{
	// do nothing
}

void update_endpoint_case_five_naive(case_five_index_naive &I, const anchor_t &a_j, const ai_t j_cost)
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

ai_t compute_case_five_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs)
{
	ai_t cost = std::numeric_limits<ai_t>::max();

	const ai_t i_a = get<0>(anchors[i]);
	const ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
	const ai_t i_c = get<1>(anchors[i]);
	const ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

	// anchor j < anchor i
	for(ai_t j = i - 1; j > 0; j--) { // dummy start anchor is handled separately
		const ai_t j_a = get<0>(anchors[j]);
		const ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]);
		const ai_t j_c = get<1>(anchors[j]);
		const ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]);

		if (costs[j] < std::numeric_limits<ai_t>::max() and
		    j_b <= i_a and i_c < j_d and j_d <= i_d) // case 5
			cost = min(cost, costs[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d));
	}

	return cost;
}

struct case_five_index init_case_five(ai_t Tlength, ai_t Qlength)
{
	return case_five_index({ MinSegmentTree<ai_t>(0, Qlength) });
}

ai_t compute_case_five(const case_five_index &I, const anchor_t &a_i)
{
	const ai_t i_a = get<0>(a_i);
	const ai_t i_c = get<1>(a_i);
	const ai_t i_d = get<1>(a_i) + get<2>(a_i);
	const ai_t i_diag = i_a - i_c;

	const ai_t rec_min = I.recursive_values.query(i_c, i_d - 1);

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

void  update_endpoint_case_five (case_five_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	const ai_t j_a = get<0>(a_j);
	const ai_t j_c = get<1>(a_j);
	const ai_t j_d = get<1>(a_j) + get<2>(a_j);
	const ai_t j_diag = j_a - j_c;

	I.recursive_values.update(j_d - 1, j_cost - j_diag);
}

struct case_two_index init_case_two(const vector<anchor_t> &anchors)
{
	list<l_t> delimiting_lines = {{0, horizontal, std::numeric_limits<ai_t>::max()}, {anchors.size()-1, horizontal, std::numeric_limits<ai_t>::min()}};
	vector<list<l_t>::iterator> anchor_diag_to_list(anchors.size());
	vector<list<l_t>::iterator> anchor_hor_to_list (anchors.size());
	anchor_hor_to_list[0] = delimiting_lines.begin();
	anchor_hor_to_list[anchors.size()-1] = --delimiting_lines.end();
	return case_two_index({
		anchors,
		std::move(delimiting_lines),
		std::move(anchor_diag_to_list),
		std::move(anchor_hor_to_list),
		{},
		{{0, 0}, {get<1>(anchors.back()) + get<2>(anchors.back()), anchors.size() - 1}},
		{}
			});
}

void prune_shadowed_delimiting_lines(case_two_index &I, const ai_t sweeping_line_a)
{
	const auto bound = I.updates.upper_bound(sweeping_line_a);
	for (auto e = I.updates.begin(); e != bound; ++e) { // TODO is this linear-time complexity or linearithmic?
		const ai_t j = e->second;
		const ai_t j_diag = get<0>(I.anchors[j]) - get<1>(I.anchors[j]);

		// check if diagonal line is still in the list or if we just added a new update before bound
		if (!I.active_diagonal_lines.contains(j_diag) or
				I.active_diagonal_lines[j_diag] != j or
				sweeping_line_a < e->first)
		       	continue;

		const auto line_it = I.anchor_diag_to_list[j];
		const auto next_it = std::next(line_it);
		const auto [_j, j_line, j_val] = *line_it;
		const auto [k,  k_line, k_val] = *next_it;
		const ai_t k_d = get<1>(I.anchors[k]) + get<2>(I.anchors[k]);
		const ai_t k_diag = get<0>(I.anchors[k]) - get<1>(I.anchors[k]);
		assert(j == _j);

		// check if update is still valid for this exact diagonal/horizontal position
		if (k_line != horizontal or k_d + j_diag != e->first)
			continue;

		if (get<2>(*std::prev(line_it)) <= k_val) { // TODO is this the correct order?
			// diagonal shadows horizontal
			cerr << "DEBUG: line " << j << "diag" << " shadows " << k << "hor" << endl;
			get<2>(*line_it) = k_val;
			I.delimiting_lines.erase(next_it);
			assert(I.active_horizontal_lines[k_d] == k);
			I.active_horizontal_lines.erase(k_d);
			const auto [h, h_line, _] = *std::next(line_it);
			if (h_line == horizontal) {
				const ai_t h_d = get<1>(I.anchors[h]) + get<2>(I.anchors[h]);
				I.updates.insert({ h_d + j_diag, j });
			}
		} else {
			// horizontal shadows diagonal
			cerr << "DEBUG: line " << j << "diag" << " is shadowed by " << k << "hor" << endl;
			I.delimiting_lines.erase(line_it);
			assert(I.active_diagonal_lines[j_diag] == j);
			I.active_diagonal_lines.erase(j_diag);

			const auto [h, h_line, _] = *std::prev(next_it);
			if (h_line == diagonal) {
				const ai_t h_diag = get<0>(I.anchors[h]) - get<1>(I.anchors[h]);
				I.updates.insert({ k_d + h_diag, h });
			}
		}
	}

	// delete events
	I.updates.erase(I.updates.begin(), I.updates.upper_bound(sweeping_line_a)); // recompute bound? TODO check
}

ai_t compute_case_two(case_two_index &I, const anchor_t &anchor_i)
{
	const ai_t i_a    = get<0>(anchor_i);
	const ai_t i_c    = get<1>(anchor_i);
	const ai_t i_diag = i_a - i_c;

	prune_shadowed_delimiting_lines(I, i_a);

	// locate closest lines
	const auto hor_lb  = std::prev(I.active_horizontal_lines.upper_bound(i_c));
	const auto diag_lb = I.active_diagonal_lines.upper_bound(i_diag); // NB we exclude the diagonal itself

	const ai_t c_hor_lb = hor_lb->first;
	const ai_t c_diag_lb = ((diag_lb != I.active_diagonal_lines.end()) ? i_a - (get<0>(I.anchors[diag_lb->second]) - get<1>(I.anchors[diag_lb->second])) : std::numeric_limits<ai_t>::min());
	ai_t rec_value;
	if (c_hor_lb <= c_diag_lb) {
		rec_value = get<2>(*(I.anchor_diag_to_list[diag_lb->second]));
	} else {
		rec_value = get<2>(*(I.anchor_hor_to_list[hor_lb->second]));
	}

	if (rec_value == std::numeric_limits<ai_t>::max()) {
		return std::numeric_limits<ai_t>::max();
	} else {
		return rec_value + i_a;
	}
}

void update_startpoint_case_two(case_two_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	// do nothing
}

void  update_endpoint_case_two (case_two_index &I, ai_t j, const anchor_t &a_j, const ai_t j_cost)
{
	const ai_t j_b = get<0>(a_j) + get<2>(a_j);
	const ai_t j_d = get<1>(a_j) + get<2>(a_j);
	const ai_t j_diag = j_b - j_d;

	prune_shadowed_delimiting_lines(I, j_b);

	// value of the two lines to (potentially?) insert
	ai_t recursion_value = j_cost - j_b;

	// find position in delimiting_lines
	list<l_t>::iterator pos;
	const auto hor_lb  = std::prev(I.active_horizontal_lines.upper_bound(j_d));
	const auto diag_lb = I.active_diagonal_lines.upper_bound(j_diag - 1); // NB we include the diagonal
	const ai_t c_hor_lb = hor_lb->first;
	const ai_t c_diag_lb = ((diag_lb != I.active_diagonal_lines.end()) ? j_b - (get<0>(I.anchors[diag_lb->second]) - get<1>(I.anchors[diag_lb->second])) : std::numeric_limits<ai_t>::min());

	if (c_hor_lb <= c_diag_lb) {
		pos = I.anchor_diag_to_list[diag_lb->second];
	} else {
		pos = I.anchor_hor_to_list [hor_lb->second];
	}
	const auto [p, p_line, p_val] = *pos;
	const ai_t p_d = get<1>(I.anchors[p]) + get<2>(I.anchors[p]);
	const ai_t p_diag = get<0>(I.anchors[p]) - get<1>(I.anchors[p]);

	if (((p_line == horizontal) ? p_d : j_b - p_diag) < j_d) {
		// case 1: no line intersects with the two new ones
		const auto j_diag_pos = I.delimiting_lines.insert(std::next(pos), { j, diagonal, p_val });
		I.anchor_diag_to_list[j] = j_diag_pos;
		I.active_diagonal_lines[j_diag] = j;
		const auto j_hor_pos =  I.delimiting_lines.insert(j_diag_pos, { j, horizontal, j_cost - j_b });
		I.anchor_hor_to_list[j] = j_hor_pos;
		I.active_horizontal_lines[j_d] = j;
		if (p_line == diagonal) {
			I.updates.insert({ j_d + p_diag, p });
		}
		const auto [k, k_line, _] = *std::next(j_diag_pos);
		if (k_line == horizontal) {
			const ai_t k_d = get<1>(I.anchors[k]) + get<2>(I.anchors[k]);
			I.updates.insert({ k_d + j_diag, j });
		}
	} else {
		// case 2: some line intersects with the two new ones
		assert(((p_line == horizontal) ? p_d : j_b - p_diag) == j_d and j_cost - j_b <= p_val);
		const auto j_diag_pos = I.delimiting_lines.insert(std::next(pos), { j, diagonal, p_val });
		I.anchor_diag_to_list[j] = j_diag_pos;
		I.active_diagonal_lines[j_diag] = j;
		const auto j_hor_pos =  I.delimiting_lines.insert(j_diag_pos, { j, horizontal, j_cost - j_b });
		I.anchor_hor_to_list[j] = j_hor_pos;
		I.active_horizontal_lines[j_d] = j;

		I.delimiting_lines.erase(pos);

		const auto [k, k_line, _] = *std::next(j_diag_pos);
		if (k_line == horizontal) {
			const ai_t k_d = get<1>(I.anchors[k]) + get<2>(I.anchors[k]);
			I.updates.insert({ k_d + j_diag, j });
		}
		assert((k_line == horizontal) ? (get<1>(I.anchors[k]) + get<2>(I.anchors[k]) != j_d) : (get<0>(I.anchors[k]) - get<1>(I.anchors[k]) != j_diag));
		// TODO possible duplication?
		const auto [h, h_line, __] = *std::prev(j_hor_pos);
		if (h_line == diagonal) {
			const ai_t h_diag = get<0>(I.anchors[h]) - get<1>(I.anchors[h]);
			I.updates.insert({ j_d + h_diag, h });
		}
	}
}

export
void chainx_backtrack(
		const vector<anchor_t> &anchors,
		const vector<ai_t> &costs,
		const chaining_mode m,
		vector<anchor_t> &chain_out
) {
	const ai_t n = anchors.size();
	chain_out.clear();

	if (costs[n-1] == std::numeric_limits<ai_t>::max())
		return;

	ai_t i = n - 1;
	chain_out.push_back(anchors[n-1]);
	if (m == semiglobal) {
		const ai_t final_c = get<1>(anchors[n-1]);
		for (ai_t j = i - 1; j > 0; j--) {
			if (costs[j] + connect_Qgap(anchors[j], final_c) == costs[n-1]) {
				chain_out.push_back(anchors[j]);
				i = j;
				break;
			}
		}
		assert(i != n - 1);
	}

	while (i > 0) {
		bool success = false;
		for (ai_t j = i - 1; j > 0; j--) {
			if (costs[j] + connect(anchors[j], anchors[i]) == costs[i] and
					chainx_precedes(anchors[j], anchors[i])) {
				chain_out.push_back(anchors[j]);
				i = j;
				success = true;
				break;
			}
		}

		if (!success) {
			if (m == global) {
				if (costs[i] == connect(anchors[0], anchors[i])) {
					break;
				}
			} else {
				if (costs[i] == connect_Qgap(anchors[0], anchors[i])) {
					break;
				}
			}
			assert(false);
		}
	}
	chain_out.push_back(anchors[0]);
	std::reverse(chain_out.begin(), chain_out.end());
}

export
void weak_backtrack(
		const vector<anchor_t> &anchors,
		const vector<ai_t> &costs,
		const chaining_mode m,
		vector<anchor_t> &chain_out
) {
	const ai_t n = anchors.size();
	chain_out.clear();

	if (costs[n-1] == std::numeric_limits<ai_t>::max())
		return;

	ai_t i = n - 1;
	chain_out.push_back(anchors[n-1]);
	if (m == semiglobal) {
		const ai_t final_c = get<1>(anchors[n-1]);
		for (ai_t j = i - 1; j > 0; j--) {
			if (costs[j] + connect_Qgap(anchors[j], final_c) == costs[n-1]) {
				chain_out.push_back(anchors[j]);
				i = j;
				break;
			}
		}
		assert(i != n - 1);
	}

	while (i > 0) {
		bool success = false;
		for (ai_t j = i - 1; j > 0; j--) {
			if (costs[j] + connect(anchors[j], anchors[i]) == costs[i] and
					weak_precedes(anchors[j], anchors[i])) {
				chain_out.push_back(anchors[j]);
				i = j;
				success = true;
				break;
			}
		}

		if (!success) {
			if (m == global) {
				if (costs[i] == connect(anchors[0], anchors[i])) {
					break;
				}
			} else {
				if (costs[i] == connect_Qgap(anchors[0], anchors[i])) {
					break;
				}
			}
			assert(false);
		}
	}
	chain_out.push_back(anchors[0]);
	std::reverse(chain_out.begin(), chain_out.end());
}

export
ai_t compute_chain_cost(const vector<anchor_t> &chain, const chaining_mode m)
{
	const ai_t n = chain.size();
	if (n <= 1) {
		return std::numeric_limits<ai_t>::max();
	}

	ai_t cost = 0;
	if (m == global) {
		cost += connect(chain[0], chain[1]);
	} else {
		cost += connect_Qgap(chain[0], chain[1]);
	}

	for (ai_t i = 1; i < ((m == global) ? n-1 : n-2); i++) {
		cost += connect(chain[i], chain[i+1]);
	}

	if (m == semiglobal) {
		cost += connect_Qgap(chain[n - 2], chain[n - 1]);
	}

	return cost;
}

} // namespace algo

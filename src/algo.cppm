export module algo;

import <vector>;
import <tuple>;
import <algorithm>; // std::sort
import <iostream>;
import <map>;
import <set>;
import <list>;
import <cassert>;
import <numeric>; // std::iota
import <sstream>;
import <string>;
import utils;
import MinSegmentTree;

using std::vector;
using std::tuple, std::pair, std::get;
using std::max, std::min, std::abs;
using std::cerr, std::endl;
using std::map, std::multimap;
using std::set, std::multiset;
using std::list;
using std::iota;
using std::ostream;
using std::string, std::to_string;
using namespace llchain;
using utils::chaining_mode, utils::chaining_mode::global, utils::chaining_mode::semiglobal, utils::anchor_t, utils::connect, utils::connect_Qgap, utils::chainx_precedes, utils::weak_precedes, utils::chain_stats;
typedef utils::anchor_index_t ai_t;

namespace llchain::algo {

/*
 * solves colinear chaining via DP under chainx or weak precedence, returns an
 *   optimal chain via backtracking; adapted from github.com/at-cg/ChainX
 * NB: O(n^2) time
 * NB: assumes anchors contains dummies (see place_dummy_anchors)
 * NB: assumes sorted anchors (see sort_anchors)
 */
export void chainx_solve_naive(
		const vector<anchor_t> &anchors,
		const chaining_mode m,
		vector<ai_t> &costs_out,
		vector<anchor_t> &chain_out
	);
export void weak_solve_naive(
		const vector<anchor_t> &anchors,
		const chaining_mode m,
		vector<ai_t> &costs_out,
		vector<anchor_t> &chain_out
	);

/*
 * solves weak-precedence colinear chaining via DP in O(n log n) time
 * NB: assumes anchors contains dummies (see place_dummy_anchors)
 * NB: assumes sorted anchors (see sort_anchors)
 * NB: assumes there are no perfect chains between the anchors (see
 *     utils::merge_perfect_chains)
 */
export void weak_solve_loglinear(
		const vector<anchor_t> &anchors,
		const ai_t Tlength,
		const ai_t Qlength,
		const chaining_mode m,
		vector<ai_t> &costs_out
	);

/*
 * debug version of weak_solve_loglinear that checks the correctness of all
 *   partial costs
 * NB: O(n^2) time or worse
 */
export void weak_solve_loglinear_debug(
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
 * compute the global or semiglobal cost of a chain using utils::connect
 * NB: assumes anchors contains dummies (see place_dummy_anchors)
 * NB: assumes the chain is colinear (see chainx_precedes and weak_precedes)
 */
export ai_t compute_chain_cost(const vector<anchor_t> &chain, const chaining_mode m);

// start of implementation

export
void chainx_solve_naive(
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

export
void weak_solve_naive(
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
			    weak_precedes(anchors[j], anchors[i])) { // weak precedence
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
 * compute diagonal/Q-pos rank of anchors: diagonal rank of i is r if the
 *   (r+1)-th hit diagonal by any anchor is diag(i)); Q-pos ranks of anchor
 *   [a..b),[c..d) are r for Q-pos c if c is the (r+1)-th hit Q-position by any
 *   anchor (c or d), r' for Q-pos d if d is the (r+1)-th hit Q-position
 */
pair<ai_t,vector<ai_t>> compute_diagonal_ranks(const vector<anchor_t> &anchors);
tuple<ai_t,vector<ai_t>,vector<ai_t>> compute_cd_ranks(const vector<anchor_t> &anchors);

/*
 * case 2: gap-gap case, larger Q-gap
 *   uses a 1-dimensional range minimum query structure (MinSegmentTree) where
 *   the (static) positions are diagonal ranks and the values (semi-dynamic,
 *   better values only) can be negative
 */
struct case_two_index {
	const vector<ai_t>& ranks; // anchor index -> diagonal rank
	MinSegmentTree<ai_t> recursive_values; // diagonal rank -> min(C[j] - q_e^j)
};
struct case_two_index init_case_two(ai_t max_rank, const vector<ai_t>& ranks);
ai_t compute_case_two(const case_two_index &I, ai_t i, const anchor_t &a_i);
void update_startpoint_case_two(case_two_index &I, ai_t j, const anchor_t &a_j, ai_t j_cost);
void update_endpoint_case_two(case_two_index &I, ai_t j, const anchor_t &a_j, const ai_t j_cost);

/*
 * helper function that considers all recursive cases considered by case 2
 */
ai_t compute_case_two_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs);

/*
 * case 3: gap-gap case, larger T-gap
 *   requires a complex data structure that handles all boundaries obtained
 *   by propagating case 2 recursions to the right (horizontal or diagonal
 *   lines)
 */
enum line_t { horizontal, diagonal };
typedef tuple<ai_t, line_t, ai_t> l_t; // (anchor index j, hor/diag, C[k] - t_e^k)
struct case_three_index {
	const vector<anchor_t>& anchors;
	list<l_t> delimiting_lines; // top-to-bottom ordered list of lines/rec values involved in case 2
	vector<list<l_t>::iterator> anchor_diag_to_list; // [j] contains pointer to diagonal   list element of anchor[j]
	vector<list<l_t>::iterator> anchor_hor_to_list;  // [j] contains pointer to horizontal list element of anchor[j]
	multimap<ai_t,ai_t> updates; // event T-position -> anchor j whose diagonal line might cross the next line
	map<ai_t, ai_t> active_horizontal_lines; // horizontal line -> anchor j
	map<ai_t, ai_t> active_diagonal_lines;   // diagonal line -> anchor j
};
struct case_three_index init_case_three(const vector<anchor_t> &anchors);
ai_t compute_case_three(case_three_index &I, const anchor_t &anchor_i);
void update_startpoint_case_three(case_three_index &I, const anchor_t &a_j, const ai_t j_cost);
void  update_endpoint_case_three (case_three_index &I, ai_t j, const anchor_t &a_j, const ai_t j_cost);

/*
 * helper function that considers all recursive cases considered by case 3
 */
ai_t compute_case_three_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs_out);

/*
 * helper function to handle intersection events for case 2 and to remove lines
 *   that have been crossed
 */
void prune_shadowed_delimiting_lines(case_three_index &I, const ai_t sweeping_line_a);

/*
 * preprocess anchor to find the closest overlaps to their startpoints
 */
vector<ai_t> compute_closest_T_overlap(const vector<anchor_t> &anchors, ai_t max_rank, const vector<ai_t> &ranks);
vector<ai_t> compute_closest_Q_overlap(const vector<anchor_t> &anchors, ai_t max_rank, const vector<ai_t> &ranks);

export
void weak_solve_loglinear(
		const vector<anchor_t> &anchors,
		const ai_t Tlength,
		const ai_t Qlength,
		const chaining_mode m,
		vector<ai_t> &costs_out
) {
	const ai_t n = anchors.size();
	costs_out = vector<ai_t>(n, 0);

	const auto [ max_diag_rank, diag_ranks ] = compute_diagonal_ranks(anchors);
	const auto closest_overlap_T = compute_closest_T_overlap(anchors, max_diag_rank, diag_ranks);
	const auto closest_overlap_Q = compute_closest_Q_overlap(anchors, max_diag_rank, diag_ranks);

	vector<ai_t> points; // horizontal line sweep (+i means T-startpoint of i-th anchor, -i means T-endpoint)
	points.reserve(2 * n - 3);
	for (ai_t i = 1; i < n - 1; i++) // skip both dummy anchors
		points.push_back(-i);
	for (ai_t i = 1; i < ((m == global) ? n : n-1); i++) // skip starting dummy anchor (and end anchor if in semi-global mode)
		points.push_back(i);
	//std::sort(points.begin(), points.end(),
	//	[&](const ai_t i, const ai_t j) -> bool {
	//		return (((i >= 0) ? 2*get<0>(anchors[i])+1 : 2*(get<0>(anchors[-i]) + get<2>(anchors[-i]))) <
	//		        ((j >= 0) ? 2*get<0>(anchors[j])+1 : 2*(get<0>(anchors[-j]) + get<2>(anchors[-j]))))
	//		        or
	//		        (i >= 0 and j >= 0 and get<0>(anchors[i]) == get<0>(anchors[j]) and get<1>(anchors[i]) < get<1>(anchors[j]));
	//});
	std::stable_sort(points.begin(), points.end(),
		[&anchors](const ai_t i, const ai_t j) -> bool {
			return (((i >= 0) ? get<0>(anchors[i]) : get<0>(anchors[-i]) + get<2>(anchors[-i])) <
			        ((j >= 0) ? get<0>(anchors[j]) : get<0>(anchors[-j]) + get<2>(anchors[-j])));
	});

	case_two_index   I_two   = init_case_two(max_diag_rank, diag_ranks);
	case_three_index I_three = init_case_three(anchors);

	for (ai_t point : points) {
		if (point >= 0) { // startpoint
			const ai_t i = point;
			const ai_t i_a = get<0>(anchors[i]);
			const ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
			const ai_t i_c = get<1>(anchors[i]);
			const ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);
			const ai_t closest_T = closest_overlap_T[i];
			const ai_t closest_Q = closest_overlap_Q[i];

			// compute cost[i]
			ai_t cost = std::numeric_limits<ai_t>::max();
			if (m == global) {
				cost = connect(anchors[0], i_a, i_b, i_c, i_d);
			} else {
				cost = connect_Qgap(anchors[0], i_c);
			}

			cost = min(cost, compute_case_two  (I_two, i, anchors[i]));
			cost = min(cost, compute_case_three(I_three, anchors[i]));
			if (closest_T != -1) {
				assert(weak_precedes(anchors[closest_T], anchors[i]));
				cost = min(cost, costs_out[closest_T] + connect(anchors[closest_T], anchors[i]));
			}
			if (closest_Q != -1) {
				assert(weak_precedes(anchors[closest_Q], anchors[i]));
				cost = min(cost, costs_out[closest_Q] + connect(anchors[closest_Q], anchors[i]));
			}
			costs_out[i] = cost;

			update_startpoint_case_two  (I_two, i, anchors[i], costs_out[i]);
			update_startpoint_case_three(I_three, anchors[i], costs_out[i]);
		} else { // endpoint
			ai_t i = -point;
			update_endpoint_case_two  (I_two, i, anchors[i], costs_out[i]);
			update_endpoint_case_three(I_three, i, anchors[i], costs_out[i]);
		}
	}

	if (m == semiglobal) {
		ai_t final_cost = std::numeric_limits<ai_t>::max();
		const ai_t final_c = get<1>(anchors[n-1]);
		for (ai_t j = 0; j < n - 1; j++) {
			if (costs_out[j] < std::numeric_limits<ai_t>::max()) {
				const ai_t c = costs_out[j] + connect_Qgap(anchors[j], final_c);
				final_cost = min(final_cost, c);
			}
		}
		costs_out[n-1] = final_cost;
	}
}

vector<ai_t> compute_closest_Q_overlap(const vector<anchor_t> &anchors, ai_t max_rank, const vector<ai_t> &ranks)
{
	const ai_t n = anchors.size();
	vector<ai_t> closest_overlap(anchors.size(), -1);
	vector<ai_t> points_Q; // horizontal line sweep moving down (+i means Q-startpoint of i-th anchor, -i means Q-endpoint)
	points_Q.reserve(2 * n - 2);
	for (ai_t i = 1; i < n - 1; i++) // skip both dummy anchors
		points_Q.push_back(-i);
	for (ai_t i = 1; i < n - 1; i++) // skip both dummy anchors
		points_Q.push_back(i);
	//std::sort(points_Q.begin(), points_Q.end(),
	//	[&](const ai_t i, const ai_t j) -> bool {
	//		return (((i >= 0) ? 2*get<1>(anchors[i])+1 : 2*(get<1>(anchors[-i]) + get<2>(anchors[-i]))) <
	//		        ((j >= 0) ? 2*get<1>(anchors[j])+1 : 2*(get<1>(anchors[-j]) + get<2>(anchors[-j]))));
	//});
	std::stable_sort(points_Q.begin(), points_Q.end(),
		[&anchors](const ai_t i, const ai_t j) -> bool {
		return (((i >= 0) ? get<1>(anchors[i]) : get<1>(anchors[-i]) + get<2>(anchors[-i])) <
		        ((j >= 0) ? get<1>(anchors[j]) : get<1>(anchors[-j]) + get<2>(anchors[-j])));
	});

	//typedef unsigned long long uai_t;
	//ordered::range_marking::Map<uai_t,ai_t> active(max_rank);
	map<ai_t,ai_t> active; // diag rank -> anchor index
	for (const ai_t point : points_Q) {
		if (point > 0) { // startpoint
			const ai_t i = point;
			const auto split = active.lower_bound(ranks[i]);
			if (split != active.begin()) {
				const auto predecessor = std::prev(split);
				closest_overlap[i] = predecessor->second;
			}
			//auto const r = active.predecessor(ranks[i]);
			//if (r.exists) {
			//	closest_overlap[i] = r.value;
			//}

			//active.insert(static_cast<uai_t>(ranks[i]), i);
			active.insert({ ranks[i], i });
		} else { // endpoint
			const ai_t i = -point;
			//bool res = active.erase(static_cast<uai_t>(ranks[i]));
			//assert(res);
			active.erase(ranks[i]);
		}
	}
	return closest_overlap;
}

vector<ai_t> compute_closest_T_overlap(const vector<anchor_t> &anchors, ai_t max_rank, const vector<ai_t> &ranks)
{
	const ai_t n = anchors.size();
	vector<ai_t> closest_overlap(anchors.size(), -1);
	vector<ai_t> points_T; // horizontal line sweep moving down (+i means T-startpoint of i-th anchor, -i means T-endpoint)
	points_T.reserve(2 * n - 2);
	for (ai_t i = 1; i < n - 1; i++) // skip both dummy anchors
		points_T.push_back(-i);
	for (ai_t i = 1; i < n - 1; i++) // skip both dummy anchors
		points_T.push_back(i);
	//std::sort(points_T.begin(), points_T.end(),
	//	[&](const ai_t i, const ai_t j) -> bool {
	//		return (((i >= 0) ? 2*get<0>(anchors[i])+1 : 2*(get<0>(anchors[-i]) + get<2>(anchors[-i]))) <
	//		        ((j >= 0) ? 2*get<0>(anchors[j])+1 : 2*(get<0>(anchors[-j]) + get<2>(anchors[-j]))));
	//});
	std::stable_sort(points_T.begin(), points_T.end(),
		[&anchors](const ai_t i, const ai_t j) -> bool {
			return (((i >= 0) ? get<0>(anchors[i]) : get<0>(anchors[-i]) + get<2>(anchors[-i])) <
			        ((j >= 0) ? get<0>(anchors[j]) : get<0>(anchors[-j]) + get<2>(anchors[-j])));
	});

	//typedef unsigned long long uai_t;
	//ordered::range_marking::Map<uai_t,ai_t> active(max_rank);
	map<ai_t,ai_t> active; // diag rank -> anchor index
	for (const ai_t point : points_T) {
		if (point > 0) { // startpoint
			const ai_t i = point;
			const auto successor = active.upper_bound(ranks[i]);
			if (successor != active.end()) {
				closest_overlap[i] = successor->second;
			}
			//auto const r = active.successor(ranks[i]);
			//if (r.exists) {
			//	closest_overlap[i] = r.value;
			//}

			//active.insert(static_cast<uai_t>(ranks[i]), i);
			active.insert({ ranks[i], i });
		} else { // endpoint
			const ai_t i = -point;
			//active.erase(static_cast<uai_t>(ranks[i]));
			active.erase(ranks[i]);
		}
	}
	return closest_overlap;
}

export
void weak_solve_loglinear_debug(
		const vector<anchor_t> &anchors,
		const ai_t Tlength,
		const ai_t Qlength,
		const chaining_mode m,
		vector<ai_t> &costs_out,
		const vector<ai_t> &correct_costs
) {
	ai_t n = anchors.size();
	costs_out = vector<ai_t>(n, 0);

	const auto [ max_diag_rank, diag_ranks ] = compute_diagonal_ranks(anchors);
	const auto closest_overlap_T = compute_closest_T_overlap(anchors, max_diag_rank, diag_ranks);
	const auto closest_overlap_Q = compute_closest_Q_overlap(anchors, max_diag_rank, diag_ranks);

	//cerr << "DEBUG: closest overlaps are defined as follows:\n";
	//for (ai_t i = 1; i < n - 1; i++) {
	//	cerr << i << "-th anchor (" << get<0>(anchors[i]) << "," << get<1>(anchors[i]) << "," << get<2>(anchors[i]) << ") has " << closest_overlap_Q[i] << "-th anchor as closest Q overlap\n";
	//}

	vector<ai_t> points; // horizontal line sweep (+i means T-startpoint of i-th anchor, -i means T-endpoint)
	points.reserve(2 * n - 1);
	for (ai_t i = 1; i < n - 1; i++) // skip both dummy anchors
		points.push_back(-i);
	for (ai_t i = 1; i < ((m == global) ? n : n-1); i++) // skip starting dummy anchor (and end anchor if in semi-global mode)
		points.push_back(i);
	std::sort(points.begin(), points.end(),
		[&](const ai_t i, const ai_t j) -> bool {
			return (((i >= 0) ? 2*get<0>(anchors[i])+1 : 2*(get<0>(anchors[-i]) + get<2>(anchors[-i]))) <
			        ((j >= 0) ? 2*get<0>(anchors[j])+1 : 2*(get<0>(anchors[-j]) + get<2>(anchors[-j]))))
			        or
			        (i >= 0 and j >= 0 and get<0>(anchors[i]) == get<0>(anchors[j]) and get<1>(anchors[i]) < get<1>(anchors[j]));
	});

	case_two_index   I_two   = init_case_two(max_diag_rank, diag_ranks);
	case_three_index I_three = init_case_three(anchors);

	for (ai_t point : points) {
		if (point >= 0) { // startpoint
			const ai_t i = point;
			const ai_t i_a = get<0>(anchors[i]);
			const ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
			const ai_t i_c = get<1>(anchors[i]);
			const ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);
			const ai_t closest_T = closest_overlap_T[i];
			const ai_t closest_Q = closest_overlap_Q[i];

			// compute cost[i]
			ai_t cost = std::numeric_limits<ai_t>::max();
			if (m == global) {
				cost = connect(anchors[0], i_a, i_b, i_c, i_d);
			} else {
				cost = connect_Qgap(anchors[0], i_c);
			}
			assert(compute_case_two  (I_two, i, anchors[i]) == compute_case_two_debug  (anchors, i, costs_out));
			cost = min(cost, compute_case_two  (I_two, i, anchors[i]));

			assert(compute_case_three(I_three,  anchors[i]) == compute_case_three_debug(anchors, i, costs_out));
			cost = min(cost, compute_case_three(I_three, anchors[i]));

			// TODO debug cases 1a and 1b
			if (closest_T != -1) {
				assert(weak_precedes(anchors[closest_T], anchors[i]));
				cost = min(cost, costs_out[closest_T] + connect(anchors[closest_T], anchors[i]));
			}
			if (closest_Q != -1) {
				assert(weak_precedes(anchors[closest_Q], anchors[i]));
				cost = min(cost, costs_out[closest_Q] + connect(anchors[closest_Q], anchors[i]));
			}

			costs_out[i] = cost;
			assert(costs_out[i] == correct_costs[i]);

			update_startpoint_case_two  (I_two, i, anchors[i], costs_out[i]);
			update_startpoint_case_three(I_three,  anchors[i], costs_out[i]);
		} else { // endpoint
			ai_t i = -point;
			update_endpoint_case_two  (I_two,   i, anchors[i], costs_out[i]);
			update_endpoint_case_three(I_three, i, anchors[i], costs_out[i]);
		}
	}

	if (m == semiglobal) {
		ai_t final_cost = std::numeric_limits<ai_t>::max();
		const ai_t final_c = get<1>(anchors[n-1]);
		for (ai_t j = 0; j < n - 1; j++) {
			if (costs_out[j] < std::numeric_limits<ai_t>::max()) {
				const ai_t c = costs_out[j] + connect_Qgap(anchors[j], final_c);
				final_cost = min(final_cost, c);
			}
		}
		costs_out[n-1] = final_cost;
	}
}

struct case_two_index init_case_two(ai_t max_rank, const vector<ai_t> &ranks)
{
	return { ranks, MinSegmentTree<ai_t>(0, max_rank) };
}

ai_t compute_case_two(const case_two_index &I, ai_t i, const anchor_t &a_i)
{
	const ai_t i_c = get<1>(a_i);
	const ai_t rec_min = I.recursive_values.query(I.ranks[i] + 1, I.recursive_values.maxquery);

	if (rec_min < std::numeric_limits<ai_t>::max()) {
		return i_c + rec_min;
	} else {
		return std::numeric_limits<ai_t>::max();
	}
}

void update_startpoint_case_two(case_two_index &I, ai_t j, const anchor_t &a_j, ai_t j_cost)
{
	// do nothing
}

void update_endpoint_case_two(case_two_index &I, ai_t j, const anchor_t &a_j, const ai_t j_cost)
{
	const ai_t j_d = get<1>(a_j) + get<2>(a_j);
	const ai_t rec_value = j_cost - j_d;
	const ai_t j_rank = I.ranks[j];

	I.recursive_values.update(j_rank, rec_value);
}

ai_t compute_case_two_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs)
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
		    j_b <= i_a and j_d <= i_c and j_a - j_c > i_a - i_c) // case 2
			cost = min(cost, costs[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d));
	}

	return cost;
}

ai_t compute_case_three_debug(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs)
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
		    j_b <= i_a and j_d <= i_c and j_a - j_c <= i_a - i_c) // case 3
			cost = min(cost, costs[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d));
	}

	return cost;
}

struct case_three_index init_case_three(const vector<anchor_t> &anchors)
{
	list<l_t> delimiting_lines = {{0, horizontal, std::numeric_limits<ai_t>::max()}, {anchors.size()-1, horizontal, std::numeric_limits<ai_t>::min()}};
	vector<list<l_t>::iterator> anchor_diag_to_list(anchors.size());
	vector<list<l_t>::iterator> anchor_hor_to_list (anchors.size());
	anchor_hor_to_list[0] = delimiting_lines.begin();
	anchor_hor_to_list[anchors.size()-1] = --delimiting_lines.end();
	return case_three_index({
		anchors,
		std::move(delimiting_lines),
		std::move(anchor_diag_to_list),
		std::move(anchor_hor_to_list),
		{},
		{{0, 0}, {get<1>(anchors.back()) + get<2>(anchors.back()), anchors.size() - 1}},
		{}
			});
}

void prune_shadowed_delimiting_lines(case_three_index &I, const ai_t sweeping_line_a)
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
		assert(j == _j);

		// check if update is still valid for this exact diagonal/horizontal position
		if (k_line != horizontal or k_d + j_diag != e->first)
			continue;

		if (get<2>(*std::prev(line_it)) <= k_val) { // TODO is this the correct order?
			// diagonal shadows horizontal
			//cerr << "DEBUG: line " << j << "diag" << " shadows " << k << "hor" << endl;
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
			//cerr << "DEBUG: line " << j << "diag" << " is shadowed by " << k << "hor" << endl;
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

ai_t compute_case_three(case_three_index &I, const anchor_t &anchor_i)
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

void update_startpoint_case_three(case_three_index &I, const anchor_t &a_j, const ai_t j_cost)
{
	// do nothing
}

void  update_endpoint_case_three(case_three_index &I, ai_t j, const anchor_t &a_j, const ai_t j_cost)
{
	const ai_t j_b = get<0>(a_j) + get<2>(a_j);
	const ai_t j_d = get<1>(a_j) + get<2>(a_j);
	const ai_t j_diag = j_b - j_d;

	prune_shadowed_delimiting_lines(I, j_b);

	// value of the two lines to (potentially?) insert
	const ai_t j_rec_value = j_cost - j_b;

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
		const auto j_hor_pos =  I.delimiting_lines.insert(j_diag_pos, { j, horizontal, j_rec_value });
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
		assert(((p_line == horizontal) ? p_d : j_b - p_diag) == j_d and j_rec_value <= p_val);
		//if (!(((p_line == horizontal) ? p_d : j_b - p_diag) == j_d and j_rec_value <= p_val)) {
		//	cerr << "WARNING: a rare edge case has been found involving case 2 and anchors " << p << ":[" << get<0>(I.anchors[p]) << ".." << get<0>(I.anchors[p]) + get<2>(I.anchors[p]) << "),[" << get<1>(I.anchors[p]) << ".." << get<1>(I.anchors[p]) + get<2>(I.anchors[p]) << ") and " << j << ":["  << get<0>(a_j) << ".." << get<0>(a_j) + get<2>(a_j) << "),[" << get<1>(a_j) << ".." << get<1>(a_j) + get<2>(a_j) << ")" << endl;
		//}
		const auto j_diag_pos = I.delimiting_lines.insert(std::next(pos), { j, diagonal, p_val });
		I.anchor_diag_to_list[j] = j_diag_pos;
		I.active_diagonal_lines[j_diag] = j;
		const auto j_hor_pos =  I.delimiting_lines.insert(j_diag_pos, { j, horizontal, j_rec_value });
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
	if (m == global and costs[n-1] == costs[0] + connect(anchors[0], anchors[n-1]))
		return;
	if (m == semiglobal and costs[n-1] == costs[0] + connect_Qgap(anchors[0], anchors[n-1]))
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
	if (m == global and costs[n-1] == costs[0] + connect(anchors[0], anchors[n-1]))
		return;
	if (m == semiglobal and costs[n-1] == costs[0] + connect_Qgap(anchors[0], anchors[n-1]))
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

pair<ai_t,vector<ai_t>> compute_diagonal_ranks(const vector<anchor_t> &anchors)
{
	const ai_t n = anchors.size();
	if (n == 0) return {0, vector<ai_t>()};

	vector<ai_t> anchors_diag(n);
	iota(anchors_diag.begin(), anchors_diag.end(), 0);
	std::sort(anchors_diag.begin(), anchors_diag.end(),
		[&anchors](ai_t i, ai_t j) -> bool {
			return get<0>(anchors[i]) - get<1>(anchors[i]) <
			       get<0>(anchors[j]) - get<1>(anchors[j]);
		});

	vector<ai_t> ranks(n);
	ranks[0] = 0;
	ai_t last_rank = 0, last_diag = get<0>(anchors[anchors_diag[0]]) - get<1>(anchors[anchors_diag[0]]);
	for (ai_t i = 1; i < n; i++) {
		const ai_t i_diag = get<0>(anchors[anchors_diag[i]]) - get<1>(anchors[anchors_diag[i]]);
		if (i_diag == last_diag) {
			ranks[anchors_diag[i]] = last_rank;
		} else {
			ranks[anchors_diag[i]] = ++last_rank;
			last_diag = i_diag;
		}
	}
	return { last_rank, ranks };
}

tuple<ai_t,vector<ai_t>,vector<ai_t>> compute_cd_ranks(const vector<anchor_t> &anchors)
{
	const ai_t n = anchors.size();
	if (n == 0) return {0, vector<ai_t>(), vector<ai_t>()};

	vector<ai_t> anchors_cd(2 * n);
	iota(anchors_cd.begin(), anchors_cd.begin() + n, 0);
	iota(anchors_cd.begin() + n, anchors_cd.end(), -(n-1));
	std::sort(anchors_cd.begin(), anchors_cd.end(),
		[&anchors](ai_t i, ai_t j) -> bool {
			return ((i >= 0) ? get<1>(anchors[i]) : get<1>(anchors[-i]) + get<2>(anchors[-i])) <
			       ((j >= 0) ? get<1>(anchors[j]) : get<1>(anchors[-j]) + get<2>(anchors[-j]));
		});

	vector<ai_t> c_ranks(n), d_ranks(n);
	c_ranks[0] = 0;
	d_ranks[0] = 0;
	ai_t last_rank = 0, last_cd = ((anchors_cd[0] >= 0) ? get<1>(anchors[anchors_cd[0]]) : get<1>(anchors[-anchors_cd[0]]) + get<2>(anchors[-anchors_cd[0]]));
	for (vector<ai_t>::size_type i = 1; i < anchors_cd.size(); i++) {
		const ai_t i_cd = ((anchors_cd[i] >= 0) ? get<1>(anchors[anchors_cd[i]]) : get<1>(anchors[-anchors_cd[i]]) + get<2>(anchors[-anchors_cd[i]]));
		if (i_cd == last_cd) {
			if (anchors_cd[i] >= 0)
				c_ranks[anchors_cd[i]] = last_rank;
			else
				d_ranks[-anchors_cd[i]] = last_rank;
		} else {
			last_rank += 1;
			if (anchors_cd[i] >= 0)
				c_ranks[anchors_cd[i]] = last_rank;
			else
				d_ranks[-anchors_cd[i]] = last_rank;
			last_cd = i_cd;
		}
	}
	return { last_rank, c_ranks, d_ranks };
}

/*
 * export a given colinear chain to CIGAR format
 * NB: the chain is expected to respect weak prec (see utils::weak_precedes)
 * NB: the chain is assumed to start and end with dummy anchors (see utils::place_dummy_anchors)
 * NB: we assume at least one non-dummy anchor
 * NB: to output a well-formed CIGAR string, we assume the chain does NOT contain adjacent anchors with connect cost = 0
 */
export
void write_cigar(
	const string::size_type text_length,
	const string::size_type query_length,
	const vector<anchor_t> &chain,
	const chaining_mode m,
	ostream &out
) {
	assert(chain.size() > 1);
	assert(chain.front() == anchor_t({-1, -1, 1}));
	assert(chain.back() == anchor_t({text_length, query_length, 1}));
	const ai_t n = chain.size();

	if (m == global) {
		const ai_t Tgap = get<0>(chain[1]);
		const ai_t Qgap = get<1>(chain[1]);
		if (min(Tgap, Qgap) > 0) {
			out << min(Tgap, Qgap) << "M";
		}
		if (Tgap > Qgap) {
			out << Tgap - Qgap <<  "D";
		}
		if (Tgap < Qgap) {
			out << Qgap - Tgap <<  "I";
		}
	} else { // semiglobal
		const ai_t Qgap = get<1>(chain[1]);
		if (Qgap > 0) {
			out << Qgap << "I";
		}
	}

	for (vector<anchor_t>::size_type i = 1; i < chain.size() - 2; i++) {
		assert(weak_precedes(chain[i], chain[i+1]));
		const ai_t i_length = get<2>(chain[i]);
		const ai_t i_a = get<0>(chain[i]), i_b = get<0>(chain[i]) + i_length;
		const ai_t i_c = get<1>(chain[i]), i_d = get<1>(chain[i]) + i_length;
		const ai_t j_a = get<0>(chain[i+1]);
		const ai_t j_c = get<1>(chain[i+1]);
		const ai_t Tgap = max((ai_t)0, j_a - i_b), Qgap = max((ai_t)0, j_c - i_d);
		const ai_t Tovl = max((ai_t)0, i_b - j_a), Qovl = max((ai_t)0, i_d - j_c);

		if (Tovl > 0 and Qovl > 0) {
			out << min(j_a - i_a, j_c - i_c) << "=";
			if (Tovl > Qovl) {
				out << Tovl - Qovl << "I";
			} else if (Tovl < Qovl) {
				out << Qovl - Tovl << "D";
			}
			// else perfect chain, do nothing
		} else if (Tovl > 0 and Qgap > 0) {
			out << j_a - i_a << "=";
			out << Tovl + Qgap << "I";
		} else if (Tgap > 0 and Qovl > 0) {
			out << j_c - i_c << "=";
			out << Tgap + Qovl << "D";
		} else {
			out << i_length - max(Tovl, Qovl) << "=";
			if (Tgap > 0 and Qgap > 0) {
				out << min(Tgap, Qgap) << "M";
			}
			if (Tgap < Qgap) {
				out << Qgap - Tgap << "I";
			}
			if (Tgap > Qgap) {
				out << Tgap - Qgap << "D";
			}
			if (Tovl > 0) {
				out << Tovl << "I";
			}
			if (Qovl > 0) {
				out << Qovl << "D";
			}
		}
	}

	out << get<2>(chain[n-2]) << "=";
	if (m == global) {
		const ai_t Tgap = get<0>(chain[n-1]) - (get<0>(chain[n-2]) + get<2>(chain[n-2]));
		const ai_t Qgap = get<1>(chain[n-1]) - (get<1>(chain[n-2]) + get<2>(chain[n-2]));
		if (min(Tgap, Qgap) > 0) {
			out << min(Tgap, Qgap) << "M";
		}
		if (Tgap > Qgap) {
			out << Tgap - Qgap <<  "D";
		}
		if (Tgap < Qgap) {
			out << Qgap - Tgap <<  "I";
		}
	} else { // semiglobal
		const ai_t Qgap = get<1>(chain[n-1]) - (get<1>(chain[n-2]) + get<2>(chain[n-2]));
		if (Qgap > 0) {
			out << Qgap << "I";
		}
	}
}

/*
 * see https://samtools.github.io/hts-specs/SAMv1.pdf
 */
export
void write_SAM_header(ostream &out)
{
	out << "@HD\tVN:1.6\n";
}

export
void write_SAM_text(const string &text, const string &text_id, ostream &out)
{
	out << "@SQ\tSN:" << text_id << "\tLN:" << text.length() << "\n";
}

export
void write_SAM_text(const long unsigned int text_length, const string &text_id, ostream &out)
{
	out << "@SQ\tSN:" << text_id << "\tLN:" << text_length << "\n";
}

/*
 * NB: see write_cigar
 */
export
void write_SAM_entry(
	const string::size_type text_length,
	const string &text_id,
	const string &query,
	const string &query_id,
	const bool is_query_rc,
	const bool store_SAM_sequence,
	const vector<anchor_t> &chain,
	const chaining_mode m,
	const ai_t anchored_ed,
	ostream &out
) {
	assert(chain.size() > 2);
	out << query_id << '\t'; // QNAME
	out << to_string(2 + (is_query_rc ? 16 : 0)) << "\t"; // FLAG
	out << text_id << '\t'; // RNAME
	out << ((m == global) ? 1 : get<0>(chain[1])+1) << '\t'; // POS
	out << "30\t"; // MAPQ
	write_cigar(text_length, query.length(), chain, m, out); // CIGAR
	out << '\t';
	out << "*\t"; // RNEXT
	out << "0\t"; // PNEXT
	out << "0\t"; // TLEN
	if (store_SAM_sequence) {
		string original(query);
		utils::reverse_complement(original);
		out << original << '\t'; // SEQ
	} else {
		out << '*' << '\t'; // SEQ
	}
	out << "*\t"; // QUAL
	out << "NM:i:" << anchored_ed; // NM TAG
	out << '\n';
}

/*
 * NB: see write_cigar
 */
export
void write_PAF_entry(
	const string::size_type text_length,
	const string &text_id,
	const string &query,
	const string &query_id,
	const bool is_query_rc,
	const vector<anchor_t> &chain,
	const chaining_mode m,
	const ai_t anchored_ed,
	ostream &out
) {
	assert(chain.size() > 2);
	const anchor_t first = chain[1]; // first non-dummy anchor
	const anchor_t last = chain[chain.size() - 2]; // last non-dummy anchor

	out <<         query_id; // QNAME
	out << "\t" << query.length(); // QLENGTH
	out << "\t" << 0; // QSTART (0-based)
	out << "\t" << query.length(); // QEND (0-based, open)
	out << "\t" << (is_query_rc ? "-" : "+"); // STRAND
	out << "\t" << text_id; // TNAME
	out << "\t" << text_length; // TLENGTH
	out << "\t" << ((m == semiglobal) ? to_string(get<0>(first)) : "0"); // TSTART (0-based, original strand)
	out << "\t" << ((m == semiglobal) ? to_string(get<0>(last) + get<2>(last)) : to_string(text_length)); // TEND (0-based, open, original strand)
	const auto [mut_cov, aln_len] = chain_stats(chain, m);
	out << "\t" << mut_cov; // MATCHES
	out << "\t" << aln_len; // ALIGNMENTLEN
	out << "\t" << 255; // QUAL
	out << "\t" << "cg:Z:"; write_cigar(text_length, query.length(), chain, m, out); // CIGAR
	out << "\n";
}

} // namespace algo

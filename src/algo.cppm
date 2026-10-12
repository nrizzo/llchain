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
import <iterator>;
import utils;
import MinSegmentTree;

using std::vector;
using std::tuple, std::pair, std::get, std::make_tuple;
using std::max, std::min, std::abs, std::reverse;
using std::cerr, std::endl;
using std::map, std::multimap;
using std::set, std::multiset;
using std::list;
using std::iota, std::numeric_limits;
using std::ostream;
using std::string, std::to_string;
using std::prev, std::next;
using namespace llchain;
using utils::chaining_mode, utils::chaining_mode::global, utils::chaining_mode::semiglobal, utils::chaining_mode::global_skipcost, utils::chaining_mode::semiglobal_skipcost, utils::precedence, utils::precedence::chainx_p, utils::precedence::weak_p, utils::anchor_t, utils::connect, utils::connect_Qgap, utils::chainx_precedes, utils::weak_precedes, utils::chain_stats;
typedef utils::anchor_index_t ai_t;

namespace llchain::algo {

/*
 * solves colinear chaining with L∞ gap costs and Δdiag overlap costs under
 *   the given precedence notion (see utils::chainx_precedes and
 *   utils::weak_precedes), returns the partial DP values (costs_out) and an
 *   optimal chain (chain_out) via backtracking;
 *   adapted from github.com/at-cg/ChainX
 * NB: O(n^2) time
 * NB: assumes anchors contains dummies (see place_dummy_anchors)
 * NB: assumes sorted anchors (see sort_anchors)
 */
export void solve_naive(
		const vector<anchor_t> &anchors,
		const chaining_mode m,
		const precedence p,
		vector<ai_t> &costs_out,
		vector<anchor_t> &chain_out
	);

/*
 * O(n log n)-time implementation of Algorithm 1
 *   (https://doi.org/10.4230/LIPIcs.WABI.2026.13)
 *   solves colinear chaining with L∞ gap costs and Δdiag overlap costs under
 *   weak precedence and returns the partial DP values (costs_out)
 * NB: assumes anchors contains dummies (see place_dummy_anchors)
 * NB: assumes sorted anchors (see sort_anchors)
 * NB: assumes there are no perfect chains between the anchors (see
 *     utils::merge_perfect_chains)
 */
export const bool solve_debug = true;
export template<bool debug> void weak_solve_loglinear(
		const vector<anchor_t> &anchors,
		const ai_t Tlength,
		const ai_t Qlength,
		const chaining_mode m,
		vector<ai_t> &costs_out
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

/*
 * export a given colinear chain to CIGAR format
 * NB: the chain is expected to respect weak prec (see utils::weak_precedes)
 * NB: the chain is assumed to start and end with dummy anchors (see utils::place_dummy_anchors)
 * NB: we assume at least one non-dummy anchor
 * NB: to output a well-formed CIGAR string, we assume the chain does NOT contain adjacent anchors with connect cost = 0
 */
export void write_cigar(
	const string::size_type text_length,
	const string::size_type query_length,
	const vector<anchor_t> &chain,
	const chaining_mode m,
	ostream &out
);

/*
 * see https://samtools.github.io/hts-specs/SAMv1.pdf
 */
export void write_SAM_header(ostream &out);
export void write_SAM_text(const string &text, const string &text_id, ostream &out);
export void write_SAM_text(const long unsigned int text_length, const string &text_id, ostream &out);

/*
 * NB: see write_cigar
 */
export void write_SAM_entry(
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
);

/*
 * NB: see write_cigar
 */
export void write_PAF_entry(
	const string::size_type text_length,
	const string &text_id,
	const string &query,
	const string &query_id,
	const bool is_query_rc,
	const vector<anchor_t> &chain,
	const chaining_mode m,
	const ai_t anchored_ed,
	ostream &out
);

// start of implementation

export
void solve_naive(
		const vector<anchor_t> &anchors,
		const chaining_mode m,
		const precedence p,
		vector<ai_t> &costs_out,
		vector<anchor_t> &chain_out
) {
	const ai_t n = anchors.size();
	costs_out = vector<ai_t>(n); // partial chainx-prec DP costs
	vector<ai_t> backtracks(n, -1); // index of an optimal previous anchor
	chain_out.clear();

	costs_out[0] = 0; // first dummy anchor
	for (ai_t i = 1; i < ((m == global) ? n : n-1); i++) {
		ai_t i_cost = numeric_limits<ai_t>::max();
		ai_t backtrack = -1;

		const ai_t i_a = get<0>(anchors[i]);
		const ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
		const ai_t i_c = get<1>(anchors[i]);
		const ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

		if (m == global) {
			i_cost = connect(anchors[0], i_a, i_b, i_c, i_d);
		} else if (m == semiglobal) {
			i_cost = connect_Qgap(anchors[0], i_c);
		} else if (m == global_skipcost) {
			i_cost = connect(anchors[0], i_a, i_b, i_c, i_d) / 2;
		} else {
			assert(m == semiglobal_skipcost);
			i_cost = connect_Qgap(anchors[0], i_c) / 2;
		}
		backtrack = 0;

		for (ai_t j = i - 1; j > 0; j--) { // anchor j < anchor i, first dummy anchor is treated separately
			const ai_t j_a = get<0>(anchors[j]);
			const ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]);
			const ai_t j_c = get<1>(anchors[j]);
			const ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]);

			if (costs_out[j] < numeric_limits<ai_t>::max() and
			    (p == chainx_p ? chainx_precedes(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d)
			                 : weak_precedes(j_a, j_c, i_a, i_c))) {
				const ai_t c = costs_out[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d);
				if (c < i_cost) backtrack = j;
				i_cost = min(i_cost, c);
			}
		}
		//save optimal cost at offset i
		costs_out[i] = i_cost;
		backtracks[i] = backtrack;
	}
	if (m != global) {
		ai_t final_cost = numeric_limits<ai_t>::max();
		ai_t backtrack = -1;
		for (ai_t j = 1; j < n - 1; j++) {
			if (costs_out[j] < numeric_limits<ai_t>::max()) {
				ai_t c;
				if (m == semiglobal) {
					const ai_t final_c = get<1>(anchors[n-1]);
					c = costs_out[j] + connect_Qgap(anchors[j], final_c);
				} else if (m == global_skipcost) {
					c = costs_out[j] + connect(anchors[j], anchors[n-1]) / 2;
				} else {
					assert(m == semiglobal_skipcost);
					const ai_t final_c = get<1>(anchors[n-1]);
					c = costs_out[j] + connect_Qgap(anchors[j], final_c) / 2;
				}
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
	reverse(chain_out.begin(), chain_out.end());
}

/*
 * compute diagonal/Q-pos rank of anchors: diagonal rank of i is r if the
 *   (r+1)-th hit diagonal by any anchor is diag(i)); Q-pos ranks of anchor
 *   [a..b),[c..d) are r for Q-pos c if c is the (r+1)-th hit Q-position by any
 *   anchor (c or d), r' for Q-pos d if d is the (r+1)-th hit Q-position
 */
pair<ai_t,vector<ai_t>> compute_diagonal_ranks(const vector<anchor_t> &anchors);

/*
 * case 1a: larger T-overlap
 *   uses a 1-dimensional dynamic successor data structure to find the closest
 *   anchor that T-overlaps the current one and has a bigger diagonal
 *   (i.e. above)
 *   NB: O(n) space, it uses the sweep line of T to save space
 */
struct case_onea_index {
	const vector<anchor_t>& anchors;
	const vector<ai_t>& ranks; // anchor index -> diagonal rank
	const vector<ai_t>& costs;
	map<ai_t,ai_t> overlap_T; // anchors that T-overlap the sweep line (diag rank -> anchor index)
};
struct case_onea_index init_case_onea(const vector<anchor_t>& anchors, const vector<ai_t>& ranks, const vector<ai_t>& costs);
ai_t compute_case_onea(const case_onea_index &I, ai_t i);
void update_startpoint_case_onea(case_onea_index &I, ai_t j);
void update_endpoint_case_onea(case_onea_index &I, ai_t j);

/*
 * case 1b: larger Q-overlap
 *   uses a 1-dimensional dynamic successor data structure to find the closest
 *   anchor that T-overlaps the current one and has a bigger diagonal
 *   (i.e. above)
 *   NB: it precomputes all closest anchors in O(n log n) time
 *   NB: Θ(n) space
 */
struct case_oneb_index {
	const vector<anchor_t>& anchors;
	const vector<ai_t>& costs;
	const vector<ai_t> closest_overlap_Q;
};
struct case_oneb_index init_case_oneb(const vector<anchor_t>& anchors, const ai_t max_diag_rank, const vector<ai_t>& ranks, const vector<ai_t>& costs);
ai_t compute_case_oneb(const case_oneb_index &I, const ai_t i);

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
struct case_two_index init_case_two(const ai_t max_rank, const vector<ai_t>& ranks);
ai_t compute_case_two(const case_two_index &I, ai_t i, const anchor_t &a_i);
void update_endpoint_case_two(case_two_index &I, ai_t j, const anchor_t &a_j, const ai_t j_cost);

/*
 * case 3: gap-gap case, larger T-gap
 *   requires a predecessor data structure to handle a 1D representation of the
 *   recursive cases projected onto the sweep line, and it needs to bookkeep all
 *   intersections between diagonal and horizontal lines and how their
 *   regions cover one another
 */
enum line_t { horizontal, diagonal };
typedef pair<ai_t, ai_t> index_recvalue_t; // (anchor index j, C[k] - t_e^k of the appropriate region)
typedef pair<ai_t, ai_t> cross; // (anchor index j, anchor index i) for a cross of j's diagonal with i's horizontal
struct case_three_index {
	const vector<anchor_t>& anchors;
	multimap<ai_t,cross> updates; // event T-position -> potential cross happening at this position
	map<ai_t,index_recvalue_t> active_horizontal_lines; // horizontal line -> anchor j, C[k] - t_e^k of the region below line
	map<ai_t,index_recvalue_t> active_diagonal_lines;   // diagonal line -> anchor j, C[k] - t_e^k of the region below line
};
struct case_three_index init_case_three(const vector<anchor_t> &anchors);
ai_t compute_case_three(case_three_index &I, const anchor_t &anchor_i);
void  update_endpoint_case_three (case_three_index &I, ai_t j, const anchor_t &a_j, const ai_t j_cost);

/*
 * helper function to handle intersection events for case 3
 */
void prune_shadowed_delimiting_lines(case_three_index &I, const ai_t sweeping_line_a);

/*
 * helper function to precompute the closest Q-overlapping anchor to the left
 *   of each anchor's startpoint
 */
vector<ai_t> compute_closest_Q_overlap(const vector<anchor_t> &anchors, ai_t max_rank, const vector<ai_t> &ranks);

/*
 * debug functions
 * NB: O(n) time
 */
ai_t compute_case_onea_naive (const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs);
ai_t compute_case_oneb_naive (const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs);
ai_t compute_case_two_naive  (const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs);
ai_t compute_case_three_naive(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs);

export
template <bool debug = false>
void weak_solve_loglinear(
		const vector<anchor_t> &anchors,
		const ai_t Tlength,
		const ai_t Qlength,
		const chaining_mode m,
		vector<ai_t> &costs_out
) {
	const ai_t n = anchors.size();
	costs_out = vector<ai_t>(n, 0); // n space

	const auto [ max_diag_rank, diag_ranks ] = compute_diagonal_ranks(anchors); // n + 1

	vector<ai_t> points; // points of horizontal line sweep (+i means T-startpoint of i-th anchor, -i means T-endpoint)
	points.reserve(2 * n - 3); // 2n - 3
	for (ai_t i = 1; i < n - 1; i++) // skip both dummy anchors
		points.push_back(-i);
	for (ai_t i = 1; i < ((m == global) ? n : n-1); i++) // skip starting dummy anchor (and end anchor if in non-global mode)
		points.push_back(i);

	// sort points for the line sweep
	std::stable_sort(points.begin(), points.end(),
		[&anchors](const ai_t i, const ai_t j) -> bool {
			return (((i >= 0) ? get<0>(anchors[i]) : get<0>(anchors[-i]) + get<2>(anchors[-i])) <
			        ((j >= 0) ? get<0>(anchors[j]) : get<0>(anchors[-j]) + get<2>(anchors[-j])));
	});

	case_onea_index  I_onea  = init_case_onea(anchors, diag_ranks, costs_out); // < n
	case_oneb_index  I_oneb  = init_case_oneb(anchors, max_diag_rank, diag_ranks, costs_out); // n
	case_two_index   I_two   = init_case_two(max_diag_rank, diag_ranks); // < 2n (2 x diagonal ranks)
	case_three_index I_three = init_case_three(anchors); // < 7n probably

	for (ai_t point : points) {
		if (point >= 0) { // startpoint of anchor [i_a..i_b),[i_c..i_d)
			const ai_t i = point;
			const ai_t i_a = get<0>(anchors[i]), i_b = get<0>(anchors[i]) + get<2>(anchors[i]);
			const ai_t i_c = get<1>(anchors[i]), i_d = get<1>(anchors[i]) + get<2>(anchors[i]);

			// compute cost[i]
			ai_t cost = numeric_limits<ai_t>::max();
			if (m == global) {
				cost = connect(anchors[0], i_a, i_b, i_c, i_d);
			} else if (m == semiglobal) {
				cost = connect_Qgap(anchors[0], i_c);
			} else if (m == global_skipcost) {
				cost = connect(anchors[0], i_a, i_b, i_c, i_d) / 2;
			} else {
				assert(m == semiglobal_skipcost);
				cost = connect_Qgap(anchors[0], i_c) / 2;
			}
			cost = min(cost, compute_case_onea (I_onea, i));
			if constexpr (debug) { assert(compute_case_onea (I_onea, i) == compute_case_onea_naive(anchors, i, costs_out)); };
			cost = min(cost, compute_case_oneb (I_oneb, i));
			if constexpr (debug) { assert(compute_case_oneb (I_oneb, i) == compute_case_oneb_naive(anchors, i, costs_out)); };
			cost = min(cost, compute_case_two  (I_two, i, anchors[i]));
			if constexpr (debug) { assert(compute_case_two  (I_two, i, anchors[i]) == compute_case_two_naive(anchors, i, costs_out)); };
			cost = min(cost, compute_case_three(I_three, anchors[i]));
			if constexpr (debug) { assert(compute_case_three(I_three, anchors[i]) == compute_case_three_naive(anchors, i, costs_out)); };
			costs_out[i] = cost;

			update_startpoint_case_onea (I_onea, i);
		} else { // endpoint
			ai_t i = -point;
			update_endpoint_case_onea (I_onea, i);
			update_endpoint_case_two  (I_two, i, anchors[i], costs_out[i]);
			update_endpoint_case_three(I_three, i, anchors[i], costs_out[i]);
		}
	}

	if (m != global) {
		ai_t final_cost = numeric_limits<ai_t>::max();
		for (ai_t j = 0; j < n - 1; j++) {
			if (costs_out[j] < numeric_limits<ai_t>::max()) {
				ai_t c;
				if (m == semiglobal) {
					const ai_t final_c = get<1>(anchors[n-1]);
					c = costs_out[j] + connect_Qgap(anchors[j], final_c);
				} else if (m == global_skipcost) {
					c = costs_out[j] + connect(anchors[j], anchors[n-1]) / 2;
				} else {
					assert(m == semiglobal_skipcost);
					const ai_t final_c = get<1>(anchors[n-1]);
					c = costs_out[j] + connect_Qgap(anchors[j], final_c) / 2;
				}
				final_cost = min(final_cost, c);
			}
		}
		costs_out[n-1] = final_cost;
	}
}

struct case_onea_index init_case_onea(const vector<anchor_t>& anchors, const vector<ai_t>& ranks, const vector<ai_t>& costs)
{
	return case_onea_index { anchors, ranks, costs, map<ai_t,ai_t>() };
}

ai_t compute_case_onea(const case_onea_index &I, ai_t i)
{
	const auto successor = I.overlap_T.upper_bound(I.ranks[i]);
	if (successor == I.overlap_T.end()) {
		return numeric_limits<ai_t>::max();
	} else {
		const ai_t closest_T = successor->second;
		assert(closest_T != -1 and weak_precedes(I.anchors[closest_T], I.anchors[i]));
		return I.costs[closest_T] + connect(I.anchors[closest_T], I.anchors[i]);
	}
}

void update_startpoint_case_onea(case_onea_index &I, ai_t j)
{
	I.overlap_T.insert({ I.ranks[j], j });
}

void update_endpoint_case_onea(case_onea_index &I, ai_t j)
{
	I.overlap_T.erase(I.ranks[j]);
}

struct case_oneb_index init_case_oneb(const vector<anchor_t>& anchors, const ai_t max_diag_rank, const vector<ai_t>& ranks, const vector<ai_t>& costs)
{
	return { anchors, costs, compute_closest_Q_overlap(anchors, max_diag_rank, ranks) };
}

ai_t compute_case_oneb(const case_oneb_index &I, const ai_t i)
{
	const ai_t closest_Q = I.closest_overlap_Q[i];
	if (closest_Q == -1) {
		return numeric_limits<ai_t>::max();
	} else {
		assert(weak_precedes(I.anchors[closest_Q], I.anchors[i]));
		return I.costs[closest_Q] + connect(I.anchors[closest_Q], I.anchors[i]);
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
	std::stable_sort(points_Q.begin(), points_Q.end(),
		[&anchors](const ai_t i, const ai_t j) -> bool {
		return (((i >= 0) ? get<1>(anchors[i]) : get<1>(anchors[-i]) + get<2>(anchors[-i])) <
		        ((j >= 0) ? get<1>(anchors[j]) : get<1>(anchors[-j]) + get<2>(anchors[-j])));
	});

	map<ai_t,ai_t> active; // diag rank -> anchor index
	for (const ai_t point : points_Q) {
		if (point > 0) { // startpoint
			const ai_t i = point;
			const auto split = active.lower_bound(ranks[i]);
			if (split != active.begin()) {
				const auto predecessor = prev(split);
				closest_overlap[i] = predecessor->second;
			}
			active.insert({ ranks[i], i });
		} else { // endpoint
			const ai_t i = -point;
			active.erase(ranks[i]);
		}
	}
	return closest_overlap;
}

struct case_two_index init_case_two(const ai_t max_rank, const vector<ai_t> &ranks)
{
	return { ranks, MinSegmentTree<ai_t>(0, max_rank) };
}

ai_t compute_case_two(const case_two_index &I, ai_t i, const anchor_t &a_i)
{
	const ai_t i_c = get<1>(a_i);
	const ai_t rec_min = I.recursive_values.query(I.ranks[i] + 1, I.recursive_values.maxquery);

	if (rec_min < numeric_limits<ai_t>::max()) {
		return i_c + rec_min;
	} else {
		return numeric_limits<ai_t>::max();
	}
}

void update_endpoint_case_two(case_two_index &I, ai_t j, const anchor_t &a_j, const ai_t j_cost)
{
	const ai_t j_d = get<1>(a_j) + get<2>(a_j);
	const ai_t rec_value = j_cost - j_d;
	const ai_t j_rank = I.ranks[j];

	I.recursive_values.update(j_rank, rec_value);
}

ai_t compute_case_onea_naive(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs)
{
	ai_t cost = numeric_limits<ai_t>::max();

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

		if (costs[j] < numeric_limits<ai_t>::max() and
		    j_a <= i_a and i_a < j_b and j_a - j_c > i_a - i_c) // case 1a
			cost = min(cost, costs[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d));
	}

	return cost;
}

ai_t compute_case_oneb_naive(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs)
{
	ai_t cost = numeric_limits<ai_t>::max();

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

		if (costs[j] < numeric_limits<ai_t>::max() and
		    j_c <= i_c and i_c < j_d and j_a - j_c < i_a - i_c) // case 1a
			cost = min(cost, costs[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d));
	}

	return cost;
}

ai_t compute_case_two_naive(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs)
{
	ai_t cost = numeric_limits<ai_t>::max();

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

		if (costs[j] < numeric_limits<ai_t>::max() and
		    j_b <= i_a and j_d <= i_c and j_a - j_c > i_a - i_c) // case 2
			cost = min(cost, costs[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d));
	}

	return cost;
}

ai_t compute_case_three_naive(const vector<anchor_t> &anchors, const ai_t i, const vector<ai_t> &costs)
{
	ai_t cost = numeric_limits<ai_t>::max();

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

		if (costs[j] < numeric_limits<ai_t>::max() and
		    j_b <= i_a and j_d <= i_c and j_a - j_c <= i_a - i_c) // case 3
			cost = min(cost, costs[j] + connect(j_a, j_b, j_c, j_d, i_a, i_b, i_c, i_d));
	}

	return cost;
}

struct case_three_index init_case_three(const vector<anchor_t> &anchors)
{
	// for convenience, there is an initial active horizontal line of the first dummy anchor
	return case_three_index({
		anchors,
		{}, // updates
		{{0, {0, numeric_limits<ai_t>::max()}}}, // active hor
		{} // active diag
	});
}

ai_t compute_case_three(case_three_index &I, const anchor_t &anchor_i)
{
	const ai_t i_a    = get<0>(anchor_i);
	const ai_t i_c    = get<1>(anchor_i);
	const ai_t i_diag = i_a - i_c;

	prune_shadowed_delimiting_lines(I, i_a);

	// locate closest active lines above anchor start and get their Q-coordinates (aka c value)
	const auto hor_up  = prev(I.active_horizontal_lines.upper_bound(i_c)); // horizontal line at i_c included
	const ai_t c_hor_up = hor_up->first; // the key is the Q-coordinate already

	const auto diag_up = I.active_diagonal_lines.upper_bound(i_diag); // diagonal line at i_diag excluded
	const ai_t c_diag_up = ((diag_up != I.active_diagonal_lines.end()) ? i_a - diag_up->first : numeric_limits<ai_t>::min());

	// find the correct region and its recursive value
	const ai_t rec_value = (c_hor_up <= c_diag_up) ? diag_up->second.second : hor_up->second.second;

	if (rec_value == numeric_limits<ai_t>::max())
		return numeric_limits<ai_t>::max();
	else
		return rec_value + i_a;
}

void  update_endpoint_case_three(case_three_index &I, ai_t j, const anchor_t &a_j, const ai_t j_cost)
{
	const ai_t j_b = get<0>(a_j) + get<2>(a_j);
	const ai_t j_d = get<1>(a_j) + get<2>(a_j);
	const ai_t j_diag = j_b - j_d;
	const ai_t j_rec_value = j_cost - j_b; // value of the new 1D segment to insert

	prune_shadowed_delimiting_lines(I, j_b);

	// find region that will be covered by our new 1D segment
	const auto hor_up  = prev(I.active_horizontal_lines.upper_bound(j_d));
	const auto diag_up = I.active_diagonal_lines.upper_bound(j_diag - 1); // NB we include the diagonal
	const ai_t c_hor_up = hor_up->first; // the key is the Q-coordinate already
	const ai_t c_diag_up = ((diag_up != I.active_diagonal_lines.end()) ? j_b - diag_up->first : numeric_limits<ai_t>::min());
	const auto [ p, p_line, p_val ] = [&]() {
		if (c_diag_up >= c_hor_up) {
			return make_tuple( diag_up->second.first, diagonal, diag_up->second.second );
		} else {
			return make_tuple( hor_up->second.first, horizontal, hor_up->second.second );
		}
	}();
	const ai_t p_d = get<1>(I.anchors[p]) + get<2>(I.anchors[p]);
	const ai_t p_diag = get<0>(I.anchors[p]) - get<1>(I.anchors[p]);

	// activate the new horizontal and diagonal lines, possibly replacing a line at the same hor/diag
	assert(j_rec_value <= p_val);
	const auto [ diag_new, _ ] = I.active_diagonal_lines.insert_or_assign(j_diag, std::make_pair(j, p_val));
	const auto [ hor_new, __ ] = I.active_horizontal_lines.insert_or_assign(j_d, std::make_pair(j, j_rec_value));

	// insert update if line following new diagonal is horizontal
	const auto [ n, n_line ] = [&]() {
		const auto hor_next  = next(hor_new);
		const auto diag_next = (diag_new != I.active_diagonal_lines.begin()) ? prev(diag_new) : I.active_diagonal_lines.end();

		if (hor_next == I.active_horizontal_lines.end() and diag_next == I.active_diagonal_lines.end()) {
			return make_tuple( (ai_t)0, diagonal ); // do nothing outcome
		}

		const ai_t c_hor_next = (hor_next != I.active_horizontal_lines.end()) ? hor_next->first : numeric_limits<ai_t>::max();
		const ai_t c_diag_next = ((diag_next != I.active_diagonal_lines.end()) ? j_b - diag_next->first : numeric_limits<ai_t>::max());

		// smaller Q-coordinate is closer, tie is won by horizontal
		if (c_hor_next <= c_diag_next) {
			return make_tuple( hor_next->second.first,  horizontal );
		} else {
			return make_tuple( diag_next->second.first, diagonal );
		}
	}();
	if (n_line == horizontal) {
		const ai_t n_d = get<1>(I.anchors[n]) + get<2>(I.anchors[n]);
		I.updates.insert({ n_d + j_diag, {j, n} });
	}

	// insert update if line preceeding new horizontal is diagonal
	if (((p_line == horizontal) ? p_d : j_b - p_diag) < j_d) { // case 1: no line intersects with the two new ones
		if (p_line == diagonal) {
			I.updates.insert({ j_d + p_diag, {p, j} });
		}
	} else { // case 2: some line intersects with the two new ones
		assert(((p_line == horizontal) ? p_d : j_b - p_diag) == j_d and j_rec_value <= p_val);
		const auto [ pp, pp_line ] = [&]() {
			const auto hor_prev = prev(hor_new); // always exists
			const ai_t c_hor_prev = hor_prev->first;

			const auto diag_prev = next(diag_new);
			const ai_t c_diag_prev = ((diag_prev != I.active_diagonal_lines.end()) ? j_b - diag_prev->first : numeric_limits<ai_t>::min());

			// larger Q-coordinate is closer, tie is won by diagonal
			if (c_hor_prev > c_diag_prev) {
				return make_tuple( hor_prev->second.first,  horizontal );
			} else {
				return make_tuple( diag_prev->second.first, diagonal );
			}
		}();
		if (pp_line == diagonal) {
			const ai_t pp_diag = get<0>(I.anchors[pp]) - get<1>(I.anchors[pp]);
			I.updates.insert({ j_d + pp_diag, {pp, j} });
		}
	}
}

void prune_shadowed_delimiting_lines(case_three_index &I, const ai_t sweeping_line_a)
{
	if (I.updates.size() == 0)
		return;

	// process the update queues in sorted order NB: the queues will be updated during this pass
	for (auto e = I.updates.begin(); e != I.updates.end() and e->first <= sweeping_line_a; ++e) {
		// check if diagonal involved in the cross is still valid
		const ai_t j = e->second.first;
		const ai_t j_diag = get<0>(I.anchors[j]) - get<1>(I.anchors[j]);
		const auto j_diag_it = I.active_diagonal_lines.find(j_diag);
		if (j_diag_it == I.active_diagonal_lines.end() or j_diag_it->second.first != j)
			continue;

		// check if horizontal involved in the cross is still valid
		const ai_t i = e->second.second;
		const ai_t i_d = get<1>(I.anchors[i]) + get<2>(I.anchors[i]);
		const auto i_hor_it = I.active_horizontal_lines.find(i_d);
		if (i_hor_it == I.active_horizontal_lines.end() or i_hor_it->second.first != i)
			continue;
		const ai_t i_val = i_hor_it->second.second;

		// find rec value of area above diagonal
		const auto [ p, p_line, p_val ] = [&]() {
			const auto hor_prev = prev(i_hor_it); // always exists
			const ai_t c_hor_prev = hor_prev->first;
			const ai_t j_hor_prev = hor_prev->second.first;
			const ai_t hor_prev_val = hor_prev->second.second;

			const auto diag_prev = next(j_diag_it);
			const ai_t c_diag_prev = (diag_prev != I.active_diagonal_lines.end()) ? e->first - diag_prev->first : numeric_limits<ai_t>::min();
			const ai_t j_diag_prev = (diag_prev != I.active_diagonal_lines.end()) ? diag_prev->second.first : 0;
			const ai_t diag_prev_val = (diag_prev != I.active_diagonal_lines.end()) ? diag_prev->second.second : numeric_limits<ai_t>::max();

			// larger Q-coordinate is closer, ties need to be double checked
			if (c_hor_prev > c_diag_prev) {
				return make_tuple(j_hor_prev, horizontal, hor_prev_val);
			} else if (c_diag_prev > c_hor_prev) {
				return make_tuple(j_diag_prev, diagonal, diag_prev_val);
			}
			// rare edge case of update on this T-coordinate!
			if (hor_prev_val <= diag_prev_val) { // this tie can be broken arbitrarily since each update triggers a scan for crosses
				return make_tuple(j_hor_prev, horizontal, hor_prev_val);
			} else {
				return make_tuple(j_diag_prev, diagonal, diag_prev_val);
			}
		}();

		if (p_val <= i_val) {
			// area above diagonal covers the area below horizontal
			//cerr << "DEBUG: line " << j << "diag" << " shadows " << i << "hor" << endl;
			j_diag_it->second.second = i_val;
			const auto [ n, n_line ] = [&]() {
				const auto hor_next  = next(i_hor_it);
				const auto diag_next = (j_diag_it != I.active_diagonal_lines.begin()) ? prev(j_diag_it) : I.active_diagonal_lines.end();

				if (hor_next == I.active_horizontal_lines.end() and diag_next == I.active_diagonal_lines.end()) {
					return make_tuple( (ai_t)0, diagonal ); // do nothing outcome
				}

				const ai_t c_hor_next = (hor_next != I.active_horizontal_lines.end()) ? hor_next->first : numeric_limits<ai_t>::max();
				const ai_t c_diag_next = (diag_next != I.active_diagonal_lines.end()) ? e->first - diag_next->first : numeric_limits<ai_t>::max();

				// smaller Q-coordinate is closer, tie is won by horizontal
				if (c_hor_next <= c_diag_next) {
					return make_tuple( hor_next->second.first,  horizontal );
				} else {
					return make_tuple( diag_next->second.first, diagonal );
				}
			}();
			if (n_line == horizontal) {
				const ai_t n_d = get<1>(I.anchors[n]) + get<2>(I.anchors[n]);
				I.updates.insert({ n_d + j_diag, {j, n}});
			}
			I.active_horizontal_lines.erase(i_d);
		} else {
			// area below horizontal covers area above diagonal
			//cerr << "DEBUG: line " << j << "diag" << " is shadowed by " << i << "hor" << endl;
			I.active_diagonal_lines.erase(j_diag);
			if (p_line == diagonal) {
				const ai_t p_diag = get<0>(I.anchors[p]) - get<1>(I.anchors[p]);
				I.updates.insert({ i_d + p_diag, {p, i} });
			}
		}
	}

	// delete events
	I.updates.erase(I.updates.begin(), I.updates.upper_bound(sweeping_line_a));
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

	if (costs[n-1] == numeric_limits<ai_t>::max())
		return;
	if (m == global and costs[n-1] == costs[0] + connect(anchors[0], anchors[n-1]))
		return;
	if (m == semiglobal and costs[n-1] == costs[0] + connect_Qgap(anchors[0], anchors[n-1]))
		return;
	if (m == global_skipcost and costs[n-1] == costs[0] + (connect(anchors[0], anchors[n-1]) / 2))
		return;
	if (m == semiglobal_skipcost and costs[n-1] == costs[0] + (connect_Qgap(anchors[0], anchors[n-1]) / 2))
		return;

	ai_t i = n - 1;
	chain_out.push_back(anchors[n-1]);
	if (m != global) {
		const ai_t final_c = get<1>(anchors[n-1]);
		for (ai_t j = i - 1; j > 0; j--) {
			if ((m == semiglobal and costs[j] + connect_Qgap(anchors[j], final_c) == costs[n-1]) or
			    (m == global_skipcost and costs[j] + (connect(anchors[j], anchors[n-1]) / 2) == costs[n-1]) or
			    (m == semiglobal_skipcost and costs[j] + (connect_Qgap(anchors[j], final_c) / 2) == costs[n-1])) {
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
				if (costs[i] == connect(anchors[0], anchors[i]))
					break;
			} else if (m == semiglobal) {
				if (costs[i] == connect_Qgap(anchors[0], anchors[i]))
					break;
			} else if (m == global_skipcost) {
				if (costs[i] == connect(anchors[0], anchors[i]) / 2)
					break;
			} else {
				assert(m == semiglobal_skipcost);
				if (costs[i] == connect_Qgap(anchors[0], anchors[i]) / 2)
					break;
			}
			assert(false);
		}
	}
	chain_out.push_back(anchors[0]);
	reverse(chain_out.begin(), chain_out.end());
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

	if (costs[n-1] == numeric_limits<ai_t>::max())
		return;
	if (m == global and costs[n-1] == costs[0] + connect(anchors[0], anchors[n-1]))
		return;
	if (m == semiglobal and costs[n-1] == costs[0] + connect_Qgap(anchors[0], anchors[n-1]))
		return;
	if (m == global_skipcost and costs[n-1] == costs[0] + (connect(anchors[0], anchors[n-1]) / 2))
		return;
	if (m == semiglobal_skipcost and costs[n-1] == costs[0] + (connect_Qgap(anchors[0], anchors[n-1]) / 2))
		return;

	ai_t i = n - 1;
	chain_out.push_back(anchors[n-1]);
	if (m != global) {
		const ai_t final_c = get<1>(anchors[n-1]);
		for (ai_t j = i - 1; j > 0; j--) {
			if ((m == semiglobal and costs[j] + connect_Qgap(anchors[j], final_c) == costs[n-1]) or
			    (m == global_skipcost and costs[j] + (connect(anchors[j], anchors[n-1]) / 2) == costs[n-1]) or
			    (m == semiglobal_skipcost and costs[j] + (connect_Qgap(anchors[j], final_c) / 2) == costs[n-1])) {
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
				if (costs[i] == connect(anchors[0], anchors[i]))
					break;
			} else if (m == semiglobal) {
				if (costs[i] == connect_Qgap(anchors[0], anchors[i]))
					break;
			} else if (m == global_skipcost) {
				if (costs[i] == connect(anchors[0], anchors[i]) / 2)
					break;
			} else {
				assert(m == semiglobal_skipcost);
				if (costs[i] == connect_Qgap(anchors[0], anchors[i]) / 2)
					break;
			}
			assert(false);
		}
	}
	chain_out.push_back(anchors[0]);
	reverse(chain_out.begin(), chain_out.end());
}

export
ai_t compute_chain_cost(const vector<anchor_t> &chain, const chaining_mode m)
{
	const ai_t n = chain.size();
	if (n <= 1) {
		return numeric_limits<ai_t>::max();
	}

	ai_t cost = 0;
	if (m == global) {
		cost += connect(chain[0], chain[1]);
	} else if (m == semiglobal) {
		cost += connect_Qgap(chain[0], chain[1]);
	} else if (m == global_skipcost) {
		cost += connect(chain[0], chain[1]) / 2;
	} else {
		assert(m == semiglobal_skipcost);
		cost += connect_Qgap(chain[0], chain[1]) / 2;
	}

	for (ai_t i = 1; i < ((m == global) ? n-1 : n-2); i++) {
		cost += connect(chain[i], chain[i+1]);
	}

	if (m == semiglobal) {
		cost += connect_Qgap(chain[n - 2], chain[n - 1]);
	} else if (m == global_skipcost) {
		cost += connect(chain[n - 2], chain[n - 1]) / 2;
	} else if (m == semiglobal_skipcost) {
		cost += connect_Qgap(chain[n - 2], chain[n - 1]) / 2;
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

	if (m == global or m == global_skipcost) {
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
	} else {
		assert(m == semiglobal or m == semiglobal_skipcost);
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
	if (m == global or m == global_skipcost) {
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
	} else {
		assert(m == semiglobal or m == semiglobal_skipcost);
		const ai_t Qgap = get<1>(chain[n-1]) - (get<1>(chain[n-2]) + get<2>(chain[n-2]));
		if (Qgap > 0) {
			out << Qgap << "I";
		}
	}
}

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
	out << ((m == global or m == global_skipcost) ? 1 : get<0>(chain[1])+1) << '\t'; // POS
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
	out << "\t" << ((m == semiglobal or m == semiglobal_skipcost) ? to_string(get<0>(first)) : "0"); // TSTART (0-based, original strand)
	out << "\t" << ((m == semiglobal or m == semiglobal_skipcost) ? to_string(get<0>(last) + get<2>(last)) : to_string(text_length)); // TEND (0-based, open, original strand)
	const auto [mut_cov, aln_len] = chain_stats(chain, m);
	out << "\t" << mut_cov; // MATCHES
	out << "\t" << aln_len; // ALIGNMENTLEN
	out << "\t" << 255; // QUAL
	out << "\t" << "cg:Z:"; write_cigar(text_length, query.length(), chain, m, out); // CIGAR
	out << "\n";
}

} // namespace algo

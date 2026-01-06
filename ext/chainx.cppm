// adaptation of github.com/algbio/ChainX to C++ module
// TODO copyright
export module chainx;

import <vector>;
import <string>;
import <tuple>;
import <algorithm>;
import <iostream>;
import <cmath>;
import <numeric>;
import <cassert>;

using std::vector;
using std::string;
using std::pair, std::tuple, std::get;
using std::min, std::max, std::abs, std::sort, std::stable_sort;
using std::cerr;
using std::floor;

namespace chainx {

export typedef long long ai_t; // anchor index type
export typedef tuple<ai_t, ai_t, ai_t> anchor_t;
export enum mode { global, semiglobal };

/**
 * compute asymmetric anchor coverage of first sequence, given the sorted anchor list
 **/
ai_t asymmetric_coverage(const vector<anchor_t> &anchors)
{
  // we assume anchors is sorted by starting position in the query
  // we assume the first and last are dummy anchors
  ai_t cov = 0;
  for (ai_t i = 1, consumed = 0; i < anchors.size() - 1; i++) {
    const ai_t qstart = std::get<0>(anchors[i]);
    const ai_t qend = std::get<0>(anchors[i]) + std::get<2>(anchors[i]) - 1;

    if (consumed - 1 >= qend)
	    continue;

    cov += qend - max(get<0>(anchors[i]), consumed) + 1;
    consumed = qend + 1;
  }
  return cov;
}

/**
 * compute anchor-restricted edit distance using strong precedence criteria
 *   optimized to run faster using engineering trick(s), comparison mode: global
 * NB: uses ChainX precedence
 * NB: assumes the anchors in input to be sorted
 **/
export
void compute_global(
		const vector<anchor_t> &anchors,
		const ai_t bound_start,
		const float ramp_up_factor,
		vector<ai_t> &costs,
		int &revisions)
{
	ai_t n = anchors.size();
	costs = vector<ai_t>(n, 0);

	ai_t bound_redit = bound_start; // distance assumed to be <= bound_start
	revisions = 0;
	// with this assumption on upper bound of distance, a gap of >bound_redit will not be allowed between adjacent anchors

	while (true) {
		ai_t inner_loop_start = 0;

		for(ai_t j = 1; j < n; j++) {
			// compute costs[i] here
			ai_t find_min_cost = std::numeric_limits<ai_t>::max();

			ai_t j_a = get<0>(anchors[j]);
			ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]) - 1;
			ai_t j_c = get<1>(anchors[j]);
			ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]) - 1;

			// anchor i < anchor j 
			while (j_a - get<0>(anchors[inner_loop_start]) - 1 > bound_redit)
				inner_loop_start++;

			for(ai_t i = j - 1; i >= inner_loop_start; i--) {
				ai_t i_a = get<0>(anchors[i]);
				ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
				ai_t i_c = get<1>(anchors[i]);
				ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

				if (costs[i] < std::numeric_limits<ai_t>::max() and 
						i_a < j_a and i_b < j_b and
						i_c < j_c and i_d < j_d) {
					ai_t gap1 = max((ai_t)0, j_a - i_b - 1);
					ai_t gap2 = max((ai_t)0, j_c - i_d - 1);
					ai_t g = max(gap1,gap2);

					ai_t overlap1 = max((ai_t)0, i_b - j_a + 1);
					ai_t overlap2 = max((ai_t)0, i_d - j_c + 1);
					ai_t o = abs(overlap1 - overlap2);

					find_min_cost = std::min(find_min_cost, costs[i] + g + o);
				}
			}
			//save optimal cost at offset j
			costs[j] = find_min_cost;
		}

		if (costs[n-1] > bound_redit) {
			bound_redit = (ai_t)((float)bound_redit * ramp_up_factor);
			revisions++;
		}
		else
			break;
	}
}

/**
 * version of compute_global that finds the optimal chain using diagonal distance
 **/
void compute_global_optimal(
		const vector<anchor_t> &anchors,
		const ai_t bound_start,
		const float ramp_up_factor,
		vector<ai_t> &costs,
		int &revisions)
{
	// anchors are sorted by starting position in first sequence (the reference)
	const ai_t n = anchors.size();
	costs = vector<ai_t>(n, 0);
	vector<ai_t> start_endpoints; // sorted (index) list of anchors startpoints and endpoints together
	start_endpoints.reserve(2*n);
	vector<ai_t> diagonal(n, 0); // diagonal (index) of each anchor
	ai_t d = -1; // number of distinct diagonals
	vector<ai_t> diagonal_bucket_value;

	ai_t bound_redit = bound_start; // distance assumed to be <= bound_start
	revisions = 0;
	// with this assumption on upper bound of distance, a gap of >bound_redit will not be allowed between adjacent anchors

	// sort anchors by start and endpoint
	for (ai_t j = 0; j < n; j++) {
		start_endpoints.push_back(j);
	}
	for (ai_t j = 0; j < n; j++) {
		start_endpoints.push_back(-j); // negative index means endpoint
	}
	stable_sort(start_endpoints.begin(), start_endpoints.end(),
			[&](const ai_t i, const ai_t j) -> bool
			{
				return (((i >= 0) ? get<0>(anchors[i]) : get<0>(anchors[-i]) + get<2>(anchors[-i]) - 1) <
				        ((j >= 0) ? get<0>(anchors[j]) : get<0>(anchors[-j]) + get<2>(anchors[-j]) - 1));
			});

	// sort anchors by diagonal to figure out each anchor's diagonal index
	{
		vector<ai_t> diagonal_order(n);
		for (ai_t j=0; j<n; j++)
			diagonal_order[j]=j;
		sort(diagonal_order.begin(), diagonal_order.end(),
				[&](const ai_t i, const ai_t j) -> bool
				{
					return (get<0>(anchors[i]) - get<1>(anchors[i])) <
					       (get<0>(anchors[j]) - get<1>(anchors[j]));
				});
		ai_t dd = -1;
		ai_t prec_diagonal = std::numeric_limits<ai_t>::max();
		for (ai_t k = 0; k < n; k++) {
			const ai_t curr_diagonal = get<0>(anchors[diagonal_order[k]]) - get<1>(anchors[diagonal_order[k]]);
			if (curr_diagonal != prec_diagonal) {
				dd += 1;
				diagonal_bucket_value.push_back(curr_diagonal);
			}
			diagonal[diagonal_order[k]] = dd;
			prec_diagonal = curr_diagonal;
		}
		d = dd + 1;
	}

	vector<ai_t> active_anchor(d, -1); // active anchor per diagonal
	while (true) {
		ai_t inner_loop_start = 0;

		// main loop
		for(ai_t j = 2; j < 2 * n; j++) {
			if (start_endpoints[j] >= 0) { // if startpoint
				const ai_t anchorj = start_endpoints[j];
				const ai_t curr_diagonal = get<0>(anchors[anchorj]) - get<1>(anchors[anchorj]);
				ai_t find_min_cost = std::numeric_limits<ai_t>::max();
				ai_t find_min_cost_gap = std::numeric_limits<ai_t>::max();
				ai_t find_min_cost_overlap = std::numeric_limits<ai_t>::max();

				// handle start of diagonal
				if (active_anchor[diagonal[anchorj]] >= 0) {
					find_min_cost = costs[active_anchor[diagonal[anchorj]]];
				}
				active_anchor[diagonal[anchorj]] = anchorj;

				ai_t j_a = get<0>(anchors[anchorj]);
				ai_t j_b = get<0>(anchors[anchorj]) + get<2>(anchors[anchorj]) - 1;
				ai_t j_c = get<1>(anchors[anchorj]);
				ai_t j_d = get<1>(anchors[anchorj]) + get<2>(anchors[anchorj]) - 1;

				// anchor i < anchor j
				while (inner_loop_start < j and
						(start_endpoints[inner_loop_start] > 0 or 
						 j_a - (get<0>(anchors[-start_endpoints[inner_loop_start]]) + get<2>(anchors[-start_endpoints[inner_loop_start]]) - 1) - 1 > bound_redit))
					inner_loop_start++;

				// process anchors with gap in first sequence
				for (ai_t i = j - 1; i >= inner_loop_start; i--) {
					if (start_endpoints[i] > 0)
						continue;
					const ai_t anchori = -start_endpoints[i];
					ai_t i_a = get<0>(anchors[anchori]);
					ai_t i_b = get<0>(anchors[anchori]) + get<2>(anchors[anchori]) - 1;
					ai_t i_c = get<1>(anchors[anchori]);
					ai_t i_d = get<1>(anchors[anchori]) + get<2>(anchors[anchori]) - 1;

					if (costs[anchori] < std::numeric_limits<ai_t>::max() and
							i_a < j_a and i_b < j_b and
							i_c < j_c and i_d < j_d) {
						ai_t gap1 = max((ai_t)0, j_a - i_b - 1);
						ai_t gap2 = max((ai_t)0, j_c - i_d - 1);
						ai_t g = max(gap1,gap2);

						ai_t overlap1 = max((ai_t)0, i_b - j_a + 1);
						ai_t overlap2 = max((ai_t)0, i_d - j_c + 1);
						ai_t o = abs(overlap1 - overlap2);

						find_min_cost_gap = min(find_min_cost_gap, costs[anchori] + g + o);
					}
				}

				// process anchors overlapping in first sequence
				for (ai_t dd = diagonal[anchorj] + 1; dd < d; dd++) {
					const ai_t diagonal_distance = abs(curr_diagonal - diagonal_bucket_value[dd]);
					if (diagonal_distance > bound_redit)
						break;
					if (active_anchor[dd] != -1) {
						if (costs[active_anchor[dd]] < std::numeric_limits<ai_t>::max()) {
							find_min_cost_overlap = min(find_min_cost_overlap, costs[active_anchor[dd]] + diagonal_distance);
						}
					}
				}
				for (ai_t dd = diagonal[anchorj] - 1; dd >= 0; dd--) {
					const ai_t diagonal_distance = abs(curr_diagonal - diagonal_bucket_value[dd]);
					if (diagonal_distance > bound_redit)
						break;
					if (active_anchor[dd] != -1) {
						if (costs[active_anchor[dd]] < std::numeric_limits<ai_t>::max()) {
							find_min_cost_overlap = min(find_min_cost_overlap, costs[active_anchor[dd]] + diagonal_distance);
						}
					}
				}

				//save optimal cost at offset j
				costs[anchorj] = min(find_min_cost_gap, find_min_cost_overlap);
				costs[anchorj] = min(costs[anchorj], find_min_cost);
			} else { // start_endpoints[j] < 0
				 // if endpoint
				const ai_t anchorj = -start_endpoints[j];
				active_anchor[diagonal[anchorj]] = -1;
			}
		}

		if (costs[n - 1] > bound_redit) {
			bound_redit = (ai_t)((float)bound_redit * ramp_up_factor);
			revisions++;
		}
		else {
			break;
		}
	}
}

/**
 * compute anchor-restricted (semi-global) edit distance using strong precedence criteria
 *   optimized to run faster using engineering trick(s)
 * NB: uses ChainX precedence
 * NB: assumes the anchors in input to be sorted
 **/
void compute_semiglobal(
		const vector<anchor_t> &anchors,
		const ai_t bound_start,
		const float ramp_up_factor,
		vector<ai_t> &costs,
		int &revisions)
{
	ai_t n = anchors.size();
	costs = vector<ai_t>(n, 0);

	ai_t bound_redit = bound_start; // distance assumed to be <= bound_start
	revisions = 0;
	// with this assumption on upper bound of distance, a gap of >bound_redit will not be allowed between adjacent anchors

	while (true)  {
		ai_t inner_loop_start = 0;

		for(ai_t j = 1; j < n; j++) {
			// compute costs[i] here
			ai_t find_min_cost = std::numeric_limits<ai_t>::max();

			ai_t j_a = get<0>(anchors[j]);
			ai_t j_b = get<0>(anchors[j]) + get<2>(anchors[j]) - 1;
			ai_t j_c = get<1>(anchors[j]);
			ai_t j_d = get<1>(anchors[j]) + get<2>(anchors[j]) - 1;

			// anchor i < anchor j 
			while (j_a - get<0>(anchors[inner_loop_start]) - 1 > bound_redit)
				inner_loop_start++;

			{
				// always consider the first dummy anchor 
				// connection to first dummy anchor is done with modified cost to allow free gaps
				ai_t i_d = get<1>(anchors[0]) + get<2>(anchors[0]) - 1;
				ai_t qry_gap = j_c - i_d - 1;
				find_min_cost = min(find_min_cost, costs[0] + qry_gap);
			}

			// process all anchors in array for the final last dummy anchor
			if (j == n - 1)
				inner_loop_start=0;

			for(ai_t i = j-1; i >= inner_loop_start; i--) {
				ai_t i_a = get<0>(anchors[i]);
				ai_t i_b = get<0>(anchors[i]) + get<2>(anchors[i]) - 1;
				ai_t i_c = get<1>(anchors[i]);
				ai_t i_d = get<1>(anchors[i]) + get<2>(anchors[i]) - 1;

				if (costs[i] < std::numeric_limits<ai_t>::max() and
						i_a < j_a and i_b < j_b and 
						i_c < j_c and i_d < j_d) {
					ai_t gap1 = max((ai_t)0, j_a - i_b - 1);
					ai_t gap2 = max((ai_t)0, j_c - i_d - 1);

					if (j == n - 1) gap1 = 0; // modified cost for the last dummy anchor to allow free gaps
					ai_t g = max(gap1, gap2);

					ai_t overlap1 = max((ai_t)0, i_b - j_a + 1);
					ai_t overlap2 = max((ai_t)0, i_d - j_c + 1);
					ai_t o = abs(overlap1 - overlap2);

					find_min_cost = min(find_min_cost, costs[i] + g + o);
				}
			}

			// save optimal cost at offset j
			costs[j] = find_min_cost;
		}

		if (costs[n - 1] > bound_redit) {
			bound_redit = (ai_t)((float)bound_redit * ramp_up_factor);
			revisions++;
		}
		else
			break;
	}
}

/**
 * version of compute_semiglobal that finds the optimal chain using diagonal
 *   distance
 **/
void compute_semiglobal_optimal(
		const vector<anchor_t> &anchors,
		const ai_t bound_start,
		const float ramp_up_factor,
		vector<ai_t> &costs,
		int &revisions)
{
	// anchors are sorted by starting position in first sequence (the reference)
	const ai_t n = anchors.size();
	costs = vector<ai_t>(n, 0);
	vector<ai_t> start_endpoints; // sorted (index) list of anchors startpoints and endpoints together
	start_endpoints.reserve(2 * n);
	vector<ai_t> diagonal(n, 0); // diagonal (index) of each anchor
	ai_t d = -1; // number of distinct diagonals
	vector<ai_t> diagonal_bucket_value;

	ai_t bound_redit = bound_start; //distance assumed to be <= bound_start
	revisions = 0;
	// with this assumption on upper bound of distance, a gap of >bound_redit will not be allowed between adjacent anchors

	// sort anchors by start and endpoint
	for (ai_t j = 0; j < n; j++) {
		start_endpoints.push_back(j);
		start_endpoints.push_back(-j); // negative index means endpoint
	}
	// TODO stable_sort here? check with compute_global_optimal
	sort(start_endpoints.begin(), start_endpoints.end(),
			[&](const ai_t i, const ai_t j) -> bool {
				return (((i >= 0) ? get<0>(anchors[i]) : get<0>(anchors[-i]) + get<2>(anchors[-i]) - 1) <
				        ((j >= 0) ? get<0>(anchors[j]) : get<0>(anchors[-j]) + get<2>(anchors[-j]) - 1));
	});

	// sort anchors by diagonal to figure out each anchor's diagonal index
	{
		vector<ai_t> diagonal_order(n);
		for (ai_t j = 0; j < n; j++)
			diagonal_order[j]=j;
		sort(diagonal_order.begin(), diagonal_order.end(),
				[&](const ai_t i, const ai_t j) -> bool {
					return (get<0>(anchors[i]) - get<1>(anchors[i])) <
					       (get<0>(anchors[j]) - get<1>(anchors[j]));
		});
		ai_t dd = -1;
		ai_t prec_diagonal = std::numeric_limits<ai_t>::max();
		for (ai_t k = 0; k < n; k++) {
			const ai_t curr_diagonal = get<0>(anchors[diagonal_order[k]]) - get<1>(anchors[diagonal_order[k]]);
			if (curr_diagonal != prec_diagonal) {
				dd += 1;
				diagonal_bucket_value.push_back(curr_diagonal);
			}
			diagonal[diagonal_order[k]] = dd;
			prec_diagonal = curr_diagonal;
		}
		d = dd + 1;
	}

	vector<ai_t> active_anchor(d, -1); // active anchor per diagonal
	while (true) {
		ai_t inner_loop_start = 0;

		// main loop
		for(ai_t j = 2; j < 2 * n; j++) {
			if (start_endpoints[j] >= 0) { // if startpoint
				const ai_t anchorj = start_endpoints[j];
				const ai_t curr_diagonal = get<0>(anchors[anchorj]) - get<1>(anchors[anchorj]);
				ai_t find_min_cost = std::numeric_limits<ai_t>::max();
				ai_t find_min_cost_gap = std::numeric_limits<ai_t>::max();
				ai_t find_min_cost_overlap = std::numeric_limits<ai_t>::max();

				// handle start of diagonal
				if (active_anchor[diagonal[anchorj]] >= 0) {
					find_min_cost = costs[active_anchor[diagonal[anchorj]]];
				}
				active_anchor[diagonal[anchorj]] = anchorj;

				ai_t j_a = get<0>(anchors[anchorj]);
				ai_t j_b = get<0>(anchors[anchorj]) + get<2>(anchors[anchorj]) - 1;
				ai_t j_c = get<1>(anchors[anchorj]);
				ai_t j_d = get<1>(anchors[anchorj]) + get<2>(anchors[anchorj]) - 1;

				// anchor i < anchor j
				while (inner_loop_start < j and
						(start_endpoints[inner_loop_start] > 0 or j_a - (get<0>(anchors[-start_endpoints[inner_loop_start]]) + get<2>(anchors[-start_endpoints[inner_loop_start]]) - 1) - 1 > bound_redit))
					inner_loop_start++;

				{
					// always consider the first dummy anchor 
					// connection to first dummy anchor is done with modified cost to allow free gaps
					ai_t i_d = get<1>(anchors[0]) + get<2>(anchors[0]) - 1;
					ai_t qry_gap = j_c - i_d - 1;
					find_min_cost = min(find_min_cost, costs[0] + qry_gap);
				}

				//process all anchors in array for the final last dummy anchor
				if (anchorj == n-1)
					inner_loop_start=0;

				// process anchors with gap in first sequence
				for(ai_t i = j - 1; i >= inner_loop_start; i--) {
					if (start_endpoints[i] > 0)
						continue;
					const ai_t anchori = -start_endpoints[i];
					ai_t i_a = get<0>(anchors[anchori]);
					ai_t i_b = get<0>(anchors[anchori]) + get<2>(anchors[anchori]) - 1;
					ai_t i_c = get<1>(anchors[anchori]);
					ai_t i_d = get<1>(anchors[anchori]) + get<2>(anchors[anchori]) - 1;

					if (costs[anchori] < std::numeric_limits<ai_t>::max() and
							i_a < j_a and i_b < j_b and 
							i_c < j_c and i_d < j_d) {
						ai_t gap1 = max((ai_t)0, j_a - i_b - 1);
						ai_t gap2 = max((ai_t)0, j_c - i_d - 1);
						if (anchorj == n - 1)
							gap1=0; //modified cost for the last dummy anchor to allow free gaps

						ai_t g = max(gap1,gap2);

						ai_t overlap1 = max((ai_t)0, i_b - j_a + 1);
						ai_t overlap2 = max((ai_t)0, i_d - j_c + 1);
						ai_t o = abs(overlap1 - overlap2);

						find_min_cost_gap = min(find_min_cost_gap, costs[anchori] + g + o);
					}
				}

				// process anchors overlapping in first sequence
				for (ai_t dd = diagonal[anchorj] + 1; dd < d; dd++) {
					const ai_t diagonal_distance = abs(curr_diagonal - diagonal_bucket_value[dd]);
					if (diagonal_distance > bound_redit)
						break;
					if (active_anchor[dd] != -1 and costs[active_anchor[dd]] < std::numeric_limits<ai_t>::max()) {
						find_min_cost_overlap = min(find_min_cost_overlap, costs[active_anchor[dd]] + diagonal_distance);
					}
				}
				for (ai_t dd = diagonal[anchorj] - 1; dd >= 0; dd--) {
					const ai_t diagonal_distance = abs(curr_diagonal - diagonal_bucket_value[dd]);
					if (diagonal_distance > bound_redit)
						break;
					if (active_anchor[dd] != -1 and costs[active_anchor[dd]] < std::numeric_limits<ai_t>::max()) {
						find_min_cost_overlap = min(find_min_cost_overlap, costs[active_anchor[dd]] + diagonal_distance);
					}
				}

				//save optimal cost at offset j
				costs[anchorj] = min(find_min_cost_gap, find_min_cost_overlap);
				costs[anchorj] = min(costs[anchorj], find_min_cost); // edge cases
			} else { // if endpoint start_endpoints[j] < 0
				const ai_t anchorj = -start_endpoints[j];
				active_anchor[diagonal[anchorj]] = -1;
			}
		}

		if (costs[n-1] > bound_redit) {
			bound_redit = (ai_t)((float)bound_redit * ramp_up_factor);
			revisions++;
		}
		else
			break;
	}
}

/**
 * run ChainX(-opt) and output chaining costs, number of main loop iterations
 * NB: assumes dummy anchors (-1,-1,0) and (|Q|,|T|,0) are in anchors (see utils::place_dummy_anchors)
 * NB: assumes the anchors have been sorted by first component (start position in T, see utils::sort_anchors)
 */
export
void chainx(
		const vector<anchor_t> &anchors,
		const ai_t qlength,
		vector<ai_t> &costs, // output ChainX-precedence DP costs
		int &revisions, // main loop iterations
		const mode m = global, // or semiglobal
		bool optimal = true, // ChainX-opt (default) or ChainX
		bool originalmagicnumbers = false)
{
	// original ChainX magic numbers
	ai_t bound_start = 100;
	float ramp_up_factor = 4;

	if (!originalmagicnumbers) {
		ai_t cov = asymmetric_coverage(anchors);
		bound_start = max((ai_t)100, (ai_t)floor(1.1*((ai_t)qlength - cov)));
		assert(bound_start >= 0);
		ramp_up_factor = 4;
	}

	if (m == global) {
		if (optimal)
			compute_global_optimal(anchors, bound_start, ramp_up_factor, costs, revisions);
		else
			compute_global(anchors, bound_start, ramp_up_factor, costs, revisions);
	} else {
		assert(m == semiglobal);
		if (optimal)
			compute_semiglobal_optimal(anchors, bound_start, ramp_up_factor, costs, revisions);
		else
			compute_semiglobal(anchors, bound_start, ramp_up_factor, costs, revisions);
	}
}

} // namespace chainx

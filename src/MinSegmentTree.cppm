export module MinSegmentTree; // this module

import <vector>;
import <algorithm>;
import <limits>;
import <cassert>;
import utils;

using std::vector;
using std::min;
using std::numeric_limits;
using q_t = llchain::utils::anchor_index_t;

/*
 * segment tree with static topology over range [minquery..maxquery] and
 *   values in T
 */
export
template<class T>
struct MinSegmentTree {
	const q_t minquery, maxquery; // extremes included
	vector<T> tree;
	T op(const T a, const T b) const {
		return min(a,b);
	}
	const T id = numeric_limits<T>::max();
	MinSegmentTree(q_t _minquery, q_t _maxquery) :
		minquery(_minquery),
		maxquery(_maxquery),
		tree(2 * (maxquery - minquery + 2), numeric_limits<T>::max())
	{}

	void update(q_t i, T val) {
		assert(minquery <= i and i <= maxquery);
		i += -minquery + 1;
		for(tree[i + (maxquery - minquery + 2)] = val, i = (i + (maxquery - minquery + 2)) / 2; i > 0; i /= 2)
			tree[i] = op(tree[2 * i], tree[2 * i + 1]);
	}

	void remove(q_t i) {
		assert(minquery <= i and i <= maxquery);
		i += -minquery + 1;
		for(tree[i + (maxquery - minquery + 2)] = id, i = (i + (maxquery - minquery + 2)) / 2; i > 0; i /= 2)
			tree[i] = op(tree[2 * i], tree[2 * i + 1]);
	}

	T query(q_t l, q_t r) const {
		if (l > r) return id; // allow empty queries
		assert(minquery <= l and l <= r and r <= maxquery);
		l += -minquery + 1; r += -minquery + 1; // [l, r]
		T lhs = T(id), rhs = T(id);
		for(l += (maxquery - minquery + 2), r += (maxquery - minquery + 2); l < r; l >>= 1, r >>= 1) {
			if(l & 1) lhs = op(lhs, tree[l++]);
			if(!(r & 1)) rhs = op(tree[r--], rhs);
		}
		return op(l == r ? op(lhs, tree[l]) : lhs, rhs);
	}
};

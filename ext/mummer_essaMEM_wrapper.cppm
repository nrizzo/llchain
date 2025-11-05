module;
#include "mummer/sparseSA.hpp" // mummer/include
#include "sparseSA.cpp"        // mummer/include/essaMEM
#include "sssort_compact.cc"   // mummer/include/essaMEM
export module mummer_essaMEM_wrapper;

import <vector>;
import <string>;
import <tuple>;
using std::vector;
using std::string;
using std::tuple;

namespace mummer_essaMEM_wrapper {

	export typedef mummer::mummer::sparseSA sparseSA;

	export sparseSA index(const string &s, int seed_min_length)
	{
		return mummer::mummer::sparseSA(mummer::mummer::sparseSA::create_auto(s.data(), s.length(), seed_min_length, true));
	}

	export void find_MEMs(const mummer::mummer::sparseSA &sa, const string &q, int seed_min_length, vector<tuple<long long, long long, long long>> &matches_out)
	{
		auto append_matches = [&](const mummer::mummer::match_t& m) { matches_out.emplace_back(m.ref, m.query, m.len); }; //0-based coordinates
		matches_out.clear();
		sa.findMEM_each(q.data(), q.length(), seed_min_length, false, append_matches);
	}

	export void find_MUMs(const mummer::mummer::sparseSA &sa, const string &q, int seed_min_length, vector<tuple<long long, long long, long long>> &matches_out)
	{
		auto append_matches = [&](const mummer::mummer::match_t& m) { matches_out.emplace_back(m.ref, m.query, m.len); }; //0-based coordinates
		matches_out.clear();
		sa.findMUM_each(q.data(), q.length(), seed_min_length, false, append_matches);
	}

} // mummer_essaMEM_wrapper

/**
 * Author: 罗穗骞, chilli
 * Date: 2019-04-11
 * License: Unknown
 * Source: Suffix array - a powerful tool for dealing with strings
 * (Chinese IOI National team training paper, 2009)
 * Description: Builds suffix array for a string.
 * \texttt{sa[i]} is the starting index of the suffix which
 * is $i$'th in the sorted suffix array.
 * The returned vector is of size $n+1$, and \texttt{sa[0] = n}.
 * The \texttt{lcp} array contains longest common prefixes for
 * neighbouring strings in the suffix array:
 * \texttt{lcp[i] = lcp(sa[i], sa[i-1])}, \texttt{lcp[0] = 0}.
 * The input string must not contain any zero bytes.
 * Time: O(n \log n)
 * Status: stress-tested
 */
#pragma once

struct SuffixArray {
	vi sa, lcp;
	SuffixArray(string& s, int lim=256) { // or basic_string<int>
		int n = sz(s) + 1, k = 0, a, b;
		vi x(all(s)+1), y(n), ws(max(n, lim)), rank(n);
		sa = lcp = y, iota(all(sa), 0);
		for (int j = 0, p = 0; p < n; j = max(1, j * 2), lim = p) {
			p = j, iota(all(y), n - j);
			rep(i,0,n) if (sa[i] >= j) y[p++] = sa[i] - j;
			fill(all(ws), 0);
			rep(i,0,n) ws[x[i]]++;
			rep(i,1,lim) ws[i] += ws[i - 1];
			for (int i = n; i--;) sa[--ws[x[y[i]]]] = y[i];
			swap(x, y), p = 1, x[sa[0]] = 0;
			rep(i,1,n) a = sa[i - 1], b = sa[i], x[b] =
				(y[a] == y[b] && y[a + j] == y[b + j]) ? p - 1 : p++;
		}
		rep(i,1,n) rank[sa[i]] = i;
		for (int i = 0, j; i < n - 1; lcp[rank[i++]] = k)
			for (k && k--, j = sa[rank[i] - 1];
					s[i + k] == s[j + k]; k++);
	}
};

vi sa_is(const vi& s, int upper) {
	int n = sz(s); if (!n) return {};
	vi sa(n); vb ls(n); /// is suffix starting at i < suffix starting at i+1
	R0F(i,n-1) ls[i] = s[i] == s[i+1] ? ls[i+1] : s[i] < s[i+1];
	/// s-type: less than next suffix -> ls[i] = 1
	/// l-type: greater than next suffix -> ls[i] = 0
	vi sum_l(upper), sum_s(upper);
	F0R(i,n) (ls[i] ? sum_l[s[i]+1] : sum_s[s[i]])++; /// note that s[i] = upper-1 -> !ls[i]
	F0R(i,upper) {
		if (i) sum_l[i] += sum_s[i-1]; /// sum_l[i] = sum_{j=0}^{i-1}(s_j+l_j)
		sum_s[i] += sum_l[i]; /// sum_s[i] = sum_{j=0}^{i-1}s_j+sum_{j=0}^{i}l_j
	}
	auto induce = [&](const vi& lms) {
		fill(all(sa),-1);
		vi buf = sum_s;
		for (int d: lms) if (d != n) sa[buf[s[d]]++] = d; /// lms is s-type, first few ...
		buf = sum_l; sa[buf[s[n-1]]++] = n-1;
		F0R(i,n) { /// do l-type in increasing order, suf[v] > suf[v+1]
			int v = sa[i]-1;
			if (v >= 0 && !ls[v]) sa[buf[s[v]]++] = v;
		}
		buf = sum_l;
		R0F(i,n) { /// do s-type in decreasing order, suf[v] < suf[v+1]
			int v = sa[i]-1;
			if (v >= 0 && ls[v]) sa[--buf[s[v]+1]] = v; /// lms is s-type, last few ...
		}
	};
	vi lms_map(n+1,-1), lms; int m = 0;
	FOR(i,1,n) if (!ls[i-1] && ls[i]) lms_map[i]=m++, lms.pb(i);
	induce(lms); // sorts LMS prefixes
	vi sorted_lms;each(v,sa)if (lms_map[v]!=-1)sorted_lms.pb(v);
	vi rec_s(m); int rec_upper = 0; // smaller subproblem
	FOR(i,1,m) { // compare two lms substrings in sorted order
		int l = sorted_lms[i-1], r = sorted_lms[i];
		int end_l = lms_map[l]+1 < m ? lms[lms_map[l]+1] : n;
		int end_r = lms_map[r]+1 < m ? lms[lms_map[r]+1] : n;
		bool same = 0; // whether lms substrings are same
		if (end_l-l == end_r-r) {
			for (;l < end_l && s[l] == s[r]; ++l,++r);
			if (l != n && s[l] == s[r]) same = 1;
		}
		rec_s[lms_map[sorted_lms[i]]] = (rec_upper += !same);
	}
	vi rec_sa = sa_is(rec_s,rec_upper+1);
	F0R(i,m) sorted_lms[i] = lms[rec_sa[i]];
	induce(sorted_lms); // sorts LMS suffixes
	return sa;
}

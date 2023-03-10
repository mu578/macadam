//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_outpabmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_zerosmxn.h>

#ifndef MC_OUTPABMXN_H
#define MC_OUTPABMXN_H

#pragma mark - mc_outpabmxn -

MC_TARGET_FUNC void mc_outpabmxnf(const int m, const int n, const int p, const int q, const int j, const int k, float * MC_TARGET_RESTRICT c, const float * a, const float * b)
{
//!# Requires c[m x n], a[m x p] and b[n x q].
//!# c=ai*bi' i.e outer product of two column-vectors.
	int l = 0, i;
	float v, e, w;
	for (; l < n; l++) {
		w = b[(q * l) + k];
		for (i = 0; i < m; i++) {
			mc_twoproductf(a[(p * i) + j], w, &v, &e);
			c[(n * i) + l] = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_outpabmxnff(const int m, const int n, const int p, const int q, const int j, const int k, double * MC_TARGET_RESTRICT c, const float * a, const float * b)
{
//!# Requires c[m x n], a[m x p] and b[n x q].
//!# c=ai*bi' i.e outer product of two column-vectors.
	int l = 0, i;
	double v, e, w;
	for (; l < n; l++) {
		w = mc_cast(double, b[(q * l) + k]);
		for (i = 0; i < m; i++) {
			mc_twoproduct(mc_cast(double, a[(p * i) + j]), w, &v, &e);
			c[(n * i) + l] = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_outpabmxnfd(const int m, const int n, const int p, const int q, const int j, const int k, double * MC_TARGET_RESTRICT c, const float * a, const double * b)
{
//!# Requires c[m x n], a[m x p] and b[n x q].
//!# c=ai*bi' i.e outer product of two column-vectors.
	int l = 0, i;
	double v, e, w;
	for (; l < n; l++) {
		w = b[(q * l) + k];
		for (i = 0; i < m; i++) {
			mc_twoproduct(mc_cast(double, a[(p * i) + j]), w, &v, &e);
			c[(n * i) + l] = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_outpabmxndf(const int m, const int n, const int p, const int q, const int j, const int k, double * MC_TARGET_RESTRICT c, const double * a, const float * b)
{
//!# Requires c[m x n], a[m x p] and b[n x q].
//!# c=ai*bi' i.e outer product of two column-vectors.
	int l = 0, i;
	double v, e, w;
	for (; l < n; l++) {
		w = mc_cast(double, b[(q * l) + k]);
		for (i = 0; i < m; i++) {
			mc_twoproduct(a[(p * i) + j], w, &v, &e);
			c[(n * i) + l] = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_outpabmxn(const int m, const int n, const int p, const int q, const int j, const int k, double * MC_TARGET_RESTRICT c, const double * a, const double * b)
{
//!# Requires c[m x n], a[m x p] and b[n x q].
//!# c=ai*bi' i.e outer product of two column-vectors.
	int l = 0, i;
	double v, e, w;
	for (; l < n; l++) {
		w = b[(q * l) + k];
		for (i = 0; i < m; i++) {
			mc_twoproduct(a[(p * i) + j], w, &v, &e);
			c[(n * i) + l] = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_outpabmxnl(const int m, const int n, const int p, const int q, const int j, const int k, long double * MC_TARGET_RESTRICT c, const long double * a, const long double * b)
{
//!# Requires c[m x n], a[m x p] and b[n x q].
//!# c=ai*bi' i.e outer product of two column-vectors.
	int l = 0, i;
	long double v, e, w;
	for (; l < n; l++) {
		w = b[(q * l) + k];
		for (i = 0; i < m; i++) {
			mc_twoproductl(a[(p * i) + j], w, &v, &e);
			c[(n * i) + l] = v + e;
		}
	}
}

#endif /* !MC_OUTPABMXN_H */

/* EOF */
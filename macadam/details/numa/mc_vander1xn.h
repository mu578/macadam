//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_vander1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_twoproduct.h>

#ifndef MC_VANDER1XN_H
#define MC_VANDER1XN_H

#pragma mark - mc_vander1xn -

MC_TARGET_FUNC void mc_vander1xnf(const int n, const int d, float * MC_TARGET_RESTRICT a, const float * x, const int f)
{
//!# Requires a[n * (d + 1)] and x[n].
//!# Returns the Vandermonde matrix A such that its columns are powers of the vector X.
//!# With the length noted n of vector X representing the m dimension and d=degree + 1
//!# the n dimension of the output matrix A.
//!# f=0: ascending powers (default).
//!# f=1: descending powers (matlab and octave).
	int i = 0, j, k;
	float w, v, e;
	for (; i < n; i++) {
		w = 1.0f;
		for (j = 0; j < d + 1; j++) {
			k    = (f == 1) ? (((d + 1) * i) + (d - j)) :  (((d + 1) * i) + j);
			a[k] = w;
			mc_twoproductf(w, x[i], &v, &e);
			w    = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_vander1xnff(const int n, const int d, double * a, const float * x, const int f)
{
//!# Requires a[n * (d + 1)] and x[n].
//!# Returns the Vandermonde matrix A such that its columns are powers of the vector X.
//!# With the length noted n of vector X representing the m dimension and d=degree + 1
//!# the n dimension of the output matrix A.
//!# f=0: ascending powers (default).
//!# f=1: descending powers (matlab and octave).
	int i = 0, j, k;
	double w, v, e;
	for (; i < n; i++) {
		w = 1.0;
		for (j = 0; j < d + 1; j++) {
			k    = (f == 1) ? (((d + 1) * i) + (d - j)) :  (((d + 1) * i) + j);
			a[k] = w;
			mc_twoproduct(w, mc_cast(double, x[i]), &v, &e);
			w    = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_vander1xn(const int n, const int d, double * MC_TARGET_RESTRICT a, const double * x, const int f)
{
//!# Requires a[n * (d + 1)] and x[n].
//!# Returns the Vandermonde matrix A such that its columns are powers of the vector X.
//!# With the length noted n of vector X representing the m dimension and d=degree + 1
//!# the n dimension of the output matrix A.
//!# f=0: ascending powers (default).
//!# f=1: descending powers (matlab and octave).
	int i = 0, j, k;
	double w, v, e;
	for (; i < n; i++) {
		w = 1.0;
		for (j = 0; j < d + 1; j++) {
			k    = (f == 1) ? (((d + 1) * i) + (d - j)) :  (((d + 1) * i) + j);
			a[k] = w;
			mc_twoproduct(w, x[i], &v, &e);
			w    = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_vander1xnl(const int n, const int d, long double * MC_TARGET_RESTRICT a, const long double * x, const int f)
{
//!# Requires a[n * (d + 1)] and x[n].
//!# Returns the Vandermonde matrix A such that its columns are powers of the vector X.
//!# With the length noted n of vector X representing the m dimension and d=degree + 1
//!# the n dimension of the output matrix A.
//!# f=0: ascending powers (default).
//!# f=1: descending powers (matlab and octave).
	int i = 0, j, k;
	long double w, v, e;
	for (; i < n; i++) {
		w = 1.0L;
		for (j = 0; j < d + 1; j++) {
			k    = (f == 1) ? (((d + 1) * i) + (d - j)) :  (((d + 1) * i) + j);
			a[k] = w;
			mc_twoproductl(w, x[i], &v, &e);
			w    = v + e;
		}
	}
}

#endif /* !MC_VANDER1XN_H */

/* EOF */
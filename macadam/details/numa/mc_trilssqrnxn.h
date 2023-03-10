//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_trilssqrnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_TRILSSQRNXN_H
#define MC_TRILSSQRNXN_H

#pragma mark - mc_trilssqrnxn -

MC_TARGET_FUNC void mc_trilssqrnxnf(const int n, const float * a, float * sumsq, float * scale, const int f)
{
//!# Computing the scaled sum of squares of the lower-triangle of A.
//!# f=0: including main diagonal.
//!# f=1: excluding main diagonal.
	int i, j = 0;
	float t;

	*scale = 0.0f;
	*sumsq = 1.0f;

	if (n > 0) {
		for (; j < ((f == 1) ? (n - 1) : n); j++) {
			for (i = (j + ((f == 1) ? 1 : 0)); i < n; i++) {
				if (0.0f != (t = mc_fabsf(a[(n * i) + j]))) {
					if (*scale < t) {
						*sumsq = 1.0f + *sumsq * mc_raise2f(*scale / t);
						*scale = t;
					} else {
						*sumsq = *sumsq + mc_raise2f(t / *scale);
					}
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_trilssqrnxnff(const int n, const float * a, double * sumsq, double * scale, const int f)
{
//!# Computing the scaled sum of squares of the lower-triangle of A.
//!# f=0: including main diagonal.
//!# f=1: excluding main diagonal.
	int i, j = 0;
	double t;

	*scale = 0.0;
	*sumsq = 1.0;

	if (n > 0) {
		for (; j < ((f == 1) ? (n - 1) : n); j++) {
			for (i = (j + ((f == 1) ? 1 : 0)); i < n; i++) {
				if (0.0 != (t = mc_fabs(mc_cast(double, a[(n * i) + j])))) {
					if (*scale < t) {
						*sumsq = 1.0 + *sumsq * mc_raise2(*scale / t);
						*scale = t;
					} else {
						*sumsq = *sumsq + mc_raise2(t / *scale);
					}
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_trilssqrnxn(const int n, const double * a, double * sumsq, double * scale, const int f)
{
//!# Computing the scaled sum of squares of the lower-triangle of A.
//!# f=0: including main diagonal.
//!# f=1: excluding main diagonal.
	int i, j = 0;
	double t;

	*scale = 0.0;
	*sumsq = 1.0;

	if (n > 0) {
		for (; j < ((f == 1) ? (n - 1) : n); j++) {
			for (i = (j + ((f == 1) ? 1 : 0)); i < n; i++) {
				if (0.0 != (t = mc_fabs(a[(n * i) + j]))) {
					if (*scale < t) {
						*sumsq = 1.0 + *sumsq * mc_raise2(*scale / t);
						*scale = t;
					} else {
						*sumsq = *sumsq + mc_raise2(t / *scale);
					}
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_trilssqrnxnl(const int n, const long double * a, long double * sumsq, long double * scale, const int f)
{
//!# Computing the scaled sum of squares of the lower-triangle of A.
//!# f=0: including main diagonal.
//!# f=1: excluding main diagonal.
	int i, j = 0;
	long double t;

	*scale = 0.0;
	*sumsq = 1.0;

	if (n > 0) {
		for (; j < ((f == 1) ? (n - 1) : n); j++) {
			for (i = (j + ((f == 1) ? 1 : 0)); i < n; i++) {
				if (0.0L != (t = mc_fabsl(a[(n * i) + j]))) {
					if (*scale < t) {
						*sumsq = 1.0L + *sumsq * mc_raise2l(*scale / t);
						*scale = t;
					} else {
						*sumsq = *sumsq + mc_raise2l(t / *scale);
					}
				}
			}
		}
	}
}

#endif /* !MC_TRILSSQRNXN_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mssqrmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/numa/mc_meanmx1.h>

#ifndef MC_MSSQRMX1_H
#define MC_MSSQRMX1_H

#pragma mark - mc_mssqrmx1 -

MC_TARGET_FUNC void mc_mssqrmx1f(const int m, const int n, const int j, const float * a, const int b, float * mean, float * sumsq, float * scale)
{
	int i;
	float t;

	*scale = 0.0f;
	*mean  = 0.0f;
	*sumsq = 1.0f;

	if (m > 0) {
		if (m == 1) {
			*mean = a[j];
		} else {
			*mean = mc_meanmx1f(m, n, j, a, b, 1);
			for (i = 0; i < m; i++) {
				if (0.0f != (t = mc_fabsf(*mean - a[(n * i) + j]))) {
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

MC_TARGET_FUNC void mc_mssqrmx1ff(const int m, const int n, const int j, const float * a, const int b, double * mean, double * sumsq, double * scale)
{
	int i;
	double t;

	*scale = 0.0;
	*mean  = 0.0;
	*sumsq = 1.0;

	if (m > 0) {
		if (m == 1) {
			*mean = a[j];
		} else {
			*mean = mc_meanmx1ff(m, n, j, a, b, 1);
			for (i = 0; i < m; i++) {
				if (0.0f != (t = mc_fabs(*mean - mc_cast(double, a[(n * i) + j])))) {
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

MC_TARGET_FUNC void mc_mssqrmx1(const int m, const int n, const int j, const double * a, const int b, double * mean, double * sumsq, double * scale)
{
	int i;
	double t;

	*scale = 0.0L;
	*mean  = 0.0L;
	*sumsq = 1.0L;

	if (m > 0) {
		if (m == 1) {
			*mean = a[j];
		} else {
			*mean = mc_meanmx1(m, n, j, a, b, 1);
			for (i = 0; i < m; i++) {
				if (0.0f != (t = mc_fabs(*mean - a[(n * i) + j]))) {
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

MC_TARGET_FUNC void mc_mssqrmx1l(const int m, const int n, const int j, const long double * a, const int b, long double * mean, long double * sumsq, long double * scale)
{
	int i;
	long double t;

	*scale = 0.0;
	*mean  = 0.0;
	*sumsq = 1.0;

	if (m > 0) {
		if (m == 1) {
			*mean = a[j];
		} else {
			*mean = mc_meanmx1l(m, n, j, a, b, 1);
			for (i = 0; i < m; i++) {
				if (0.0f != (t = mc_fabsl(*mean - a[(n * i) + j]))) {
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

#endif /* !MC_MSSQRMX1_H */

/* EOF */
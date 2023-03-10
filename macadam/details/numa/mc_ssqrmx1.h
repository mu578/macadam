//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ssqrmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_SSQRMX1_H
#define MC_SSQRMX1_H

#pragma mark - mc_ssqrmx1 -

MC_TARGET_FUNC void mc_ssqrmx1f(const int m, const int n, const int j, const float * a, float * sumsq, float * scale)
{
	int i;
	float t;

	*scale = 0.0f;
	*sumsq = 1.0f;

	if (m > 0) {
		if (m == 1) {
			*scale = mc_fabsf(a[j]);
		} else {
			for (i = 0; i < m; i++) {
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

MC_TARGET_FUNC void mc_ssqrmx1ff(const int m, const int n, const int j, const float * a, double * sumsq, double * scale)
{
	int i;
	double t;

	*scale = 0.0;
	*sumsq = 1.0;

	if (m > 0) {
		if (m == 1) {
			*scale = mc_fabs(a[j]);
		} else {
			for (i = 0; i < m; i++) {
				if (0.0f != (t = mc_fabs(mc_cast(double, a[(n * i) + j])))) {
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

MC_TARGET_FUNC void mc_ssqrmx1(const int m, const int n, const int j, const double * a, double * sumsq, double * scale)
{
	int i;
	double t;

	*scale = 0.0L;
	*sumsq = 1.0L;

	if (m > 0) {
		if (m == 1) {
			*scale = mc_fabs(a[j]);
		} else {
			for (i = 0; i < m; i++) {
				if (0.0f != (t = mc_fabs(a[(n * i) + j]))) {
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

MC_TARGET_FUNC void mc_ssqrmx1l(const int m, const int n, const int j, const long double * a, long double * sumsq, long double * scale)
{
	int i;
	long double t;

	*scale = 0.0;
	*sumsq = 1.0;

	if (m > 0) {
		if (m == 1) {
			*scale = mc_fabsl(a[j]);
		} else {
			for (i = 0; i < m; i++) {
				if (0.0f != (t = mc_fabsl(a[(n * i) + j]))) {
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

#endif /* !MC_SSQRMX1_H */

/* EOF */
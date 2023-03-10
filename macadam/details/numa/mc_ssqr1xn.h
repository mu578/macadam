//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ssqr1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_SSQR1XN_H
#define MC_SSQR1XN_H

#pragma mark - mc_ssqr1xn -

MC_TARGET_FUNC void mc_ssqr1xnf(const int n, const float * x, float * sumsq, float * scale)
{
	int i;
	float t;

	*scale = 0.0f;
	*sumsq = 1.0f;

	if (n > 0) {
		if (n == 1) {
			*scale = mc_fabsf(x[0]);
		} else {
			for (i = 0; i < n; i++) {
				if (0.0f != (t = mc_fabsf(x[i]))) {
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

MC_TARGET_FUNC void mc_ssqr1xnff(const int n, const float * x, double * sumsq, double * scale)
{
	int i;
	double t;

	*scale = 0.0;
	*sumsq = 1.0;

	if (n > 0) {
		if (n == 1) {
			*scale = mc_fabs(x[0]);
		} else {
			for (i = 0; i < n; i++) {
				if (0.0 != (t = mc_fabs(mc_cast(double, x[i])))) {
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

MC_TARGET_FUNC void mc_ssqr1xn(const int n, const double * x, double * sumsq, double * scale)
{
	int i;
	double t;

	*scale = 0.0L;
	*sumsq = 1.0L;

	if (n > 0) {
		if (n == 1) {
			*scale = mc_fabs(x[0]);
		} else {
			for (i = 0; i < n; i++) {
				if (0.0 != (t = mc_fabs(x[i]))) {
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

MC_TARGET_FUNC void mc_ssqr1xnl(const int n, const long double * x, long double * sumsq, long double * scale)
{
	int i;
	long double t;

	*scale = 0.0;
	*sumsq = 1.0;

	if (n > 0) {
		if (n == 1) {
			*scale = mc_fabsl(x[0]);
		} else {
			for (i = 0; i < n; i++) {
				if (0.0L != (t = mc_fabsl(x[i]))) {
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

#endif /* !MC_SSQR1XN_H */

/* EOF */
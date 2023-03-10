//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ssqr2x1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_SSQR2X1_H
#define MC_SSQR2X1_H

#pragma mark - mc_ssqr2x1 -

MC_TARGET_FUNC void mc_ssqr2x1f(const int n, const int j, const float * a, float * sumsq, float * scale)
{
	float t;

	*scale = 0.0f;
	*sumsq = 1.0f;

	if (0.0f != (t = mc_fabsf(a[j]))) {
		if (*scale < t) {
			*sumsq = 1.0f + *sumsq * mc_raise2f(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2f(t / *scale);
		}
	}
	if (0.0f != (t = mc_fabsf(a[n + j]))) {
		if (*scale < t) {
			*sumsq = 1.0f + *sumsq * mc_raise2f(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2f(t / *scale);
		}
	}
}

MC_TARGET_FUNC void mc_ssqr2x1ff(const int n, const int j, const float * a, double * sumsq, double * scale)
{
	double t;

	*scale = 0.0;
	*sumsq = 1.0;

	if (0.0 != (t = mc_fabs(mc_cast(double, a[j])))) {
		if (*scale < t) {
			*sumsq = 1.0 + *sumsq * mc_raise2(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2(t / *scale);
		}
	}
	if (0.0 != (t = mc_fabs(mc_cast(double, a[n + j])))) {
		if (*scale < t) {
			*sumsq = 1.0 + *sumsq * mc_raise2(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2(t / *scale);
		}
	}
}

MC_TARGET_FUNC void mc_ssqr2x1(const int n, const int j, const double * a, double * sumsq, double * scale)
{
	double t;

	*scale = 0.0;
	*sumsq = 1.0;

	if (0.0 != (t = mc_fabs(a[j]))) {
		if (*scale < t) {
			*sumsq = 1.0 + *sumsq * mc_raise2(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2(t / *scale);
		}
	}
	if (0.0 != (t = mc_fabs(a[n + j]))) {
		if (*scale < t) {
			*sumsq = 1.0 + *sumsq * mc_raise2(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2(t / *scale);
		}
	}
}

MC_TARGET_FUNC void mc_ssqr2x1l(const int n, const int j, const long double * a, long double * sumsq, long double * scale)
{
	long double t;

	*scale = 0.0L;
	*sumsq = 1.0L;

	if (0.0L != (t = mc_fabsl(a[j]))) {
		if (*scale < t) {
			*sumsq = 1.0L + *sumsq * mc_raise2l(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2l(t / *scale);
		}
	}
	if (0.0L != (t = mc_fabsl(a[n + j]))) {
		if (*scale < t) {
			*sumsq = 1.0L + *sumsq * mc_raise2l(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2l(t / *scale);
		}
	}
}

#endif /* !MC_SSQR2X1_H */

/* EOF */
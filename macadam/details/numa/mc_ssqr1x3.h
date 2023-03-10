//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ssqr1x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_SSQR1X3_H
#define MC_SSQR1X3_H

#pragma mark - mc_ssqr1x3 -

MC_TARGET_FUNC void mc_ssqr1x3f(const float x[3], float * sumsq, float * scale)
{
	float t;

	*scale = 0.0f;
	*sumsq = 1.0f;

	if (0.0f != (t = mc_fabsf(x[0]))) {
		if (*scale < t) {
			*sumsq = 1.0f + *sumsq * mc_raise2f(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2f(t / *scale);
		}
	}
	if (0.0f != (t = mc_fabsf(x[1]))) {
		if (*scale < t) {
			*sumsq = 1.0f + *sumsq * mc_raise2f(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2f(t / *scale);
		}
	}
	if (0.0f != (t = mc_fabsf(x[2]))) {
		if (*scale < t) {
			*sumsq = 1.0f + *sumsq * mc_raise2f(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2f(t / *scale);
		}
	}
}

MC_TARGET_FUNC void mc_ssqr1x3ff(const float x[3], double * sumsq, double * scale)
{
	double t;

	*scale = 0.0;
	*sumsq = 1.0;

	if (0.0 != (t = mc_fabs(mc_cast(double, x[0])))) {
		if (*scale < t) {
			*sumsq = 1.0 + *sumsq * mc_raise2(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2(t / *scale);
		}
	}
	if (0.0 != (t = mc_fabs(mc_cast(double, x[1])))) {
		if (*scale < t) {
			*sumsq = 1.0 + *sumsq * mc_raise2(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2(t / *scale);
		}
	}
	if (0.0 != (t = mc_fabs(mc_cast(double, x[2])))) {
		if (*scale < t) {
			*sumsq = 1.0 + *sumsq * mc_raise2(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2(t / *scale);
		}
	}
}

MC_TARGET_FUNC void mc_ssqr1x3(const double x[3], double * sumsq, double * scale)
{
	double t;

	*scale = 0.0L;
	*sumsq = 1.0L;

	if (0.0 != (t = mc_fabs(x[0]))) {
		if (*scale < t) {
			*sumsq = 1.0 + *sumsq * mc_raise2(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2(t / *scale);
		}
	}
	if (0.0 != (t = mc_fabs(x[1]))) {
		if (*scale < t) {
			*sumsq = 1.0 + *sumsq * mc_raise2(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2(t / *scale);
		}
	}
	if (0.0 != (t = mc_fabs(x[2]))) {
		if (*scale < t) {
			*sumsq = 1.0 + *sumsq * mc_raise2(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2(t / *scale);
		}
	}
}

MC_TARGET_FUNC void mc_ssqr1x3l(const long double x[3], long double * sumsq, long double * scale)
{
	long double t;

	*scale = 0.0;
	*sumsq = 1.0;

	if (0.0L != (t = mc_fabsl(x[0]))) {
		if (*scale < t) {
			*sumsq = 1.0L + *sumsq * mc_raise2l(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2l(t / *scale);
		}
	}
	if (0.0L != (t = mc_fabsl(x[1]))) {
		if (*scale < t) {
			*sumsq = 1.0L + *sumsq * mc_raise2l(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2l(t / *scale);
		}
	}
	if (0.0L != (t = mc_fabsl(x[2]))) {
		if (*scale < t) {
			*sumsq = 1.0L + *sumsq * mc_raise2l(*scale / t);
			*scale = t;
		} else {
			*sumsq = *sumsq + mc_raise2l(t / *scale);
		}
	}
}

#endif /* !MC_SSQR1X3_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rempio2_cw.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fisval.h>
#include <macadam/details/math/mc_floor.h>
#include <macadam/details/math/mc_itrunc64.h>

#ifndef MC_REMPIO2_CW_H
#define MC_REMPIO2_CW_H

MC_TARGET_PROC int64_t mc_rempio2_cwf(const float x, float * z)
{
//!# Cody-Waite reduction of x such as z=x - n*pi/2 where |z|<=pi/4, result = n mod 8.
	const float dp1 = +0.78538513183593750000000000000000000000E+00f;
	const float dp2 = +1.30315311253070831298828125000000000000E-05f;
	const float dp3 = +3.03855025325309630000000000000000000000E-11f;

	int64_t r = 0;
	float w;
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	const int64_t max = INT64_MAX;
#	else
	const int64_t max = mc_cast(const int64_t, 0x7FFFFFFFFFFFFFFF);
#	endif

	*z = 0.0;
	if (mc_fisvalf(x) && mc_fabsf(x) <= 1048576.0f) {
		if (x == 0.0f) {
			return 0.0f;
		}
		w = mc_floorf(x * MCK_KF(MCK_1_PI_4));
		r = mc_itrunc64f(w - 16.0f * mc_floorf(w / 16.0f));
		if (r < (max - 1)) {
			if (r & 1) {
				++r;
				w = w + 1.0f;
			}
			 r = (mc_cast(uint64_t, r) >> 1) & 7;
			*z = ((x - w * dp1) - w * dp2) - w * dp3;
		} else {
			 r = max;
			*z = MCK_NAN;
		}
	}
	return r;
}

MC_TARGET_PROC int64_t mc_rempio2_cw(const double x, double * z)
{
//!# Cody-Waite reduction of x such as z=x - n*pi/2 where |z|<=pi/4, result = n mod 8.
	const double dp1 = +7.85398125648498535156000000000000000000000E-01;
	const double dp2 = +3.77489470793079817668000000000000000000000E-08;
	const double dp3 = +2.69515142907905952645000000000000000000000E-15;

	int64_t r = 0;
	double w;
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	const int64_t max = INT64_MAX;
#	else
	const int64_t max = mc_cast(const int64_t, 0x7FFFFFFFFFFFFFFF);
#	endif

	*z = 0.0;
	if (mc_fisval(x) && mc_fabs(x) <= 2147483648.0) {
		if (x == 0.0) {
			return 0.0;
		}
		w = mc_floor(x * MCK_K(MCK_1_PI_4));
		r = mc_itrunc64(w - 16.0 * mc_floor(w / 16.0));
		if (r < (max - 1)) {
			if (r & 1) {
				++r;
				w = w + 1.0;
			}
			 r = (mc_cast(uint64_t, r) >> 1) & 7;
			*z = ((x - w * dp1) - w * dp2) - w * dp3;
		} else {
			 r = max;
			*z = MCK_NAN;
		}
	}
	return r;
}

MC_TARGET_PROC int64_t mc_rempio2_cwl(const long double x, long double * z)
{
#	if MC_TARGET_HAVE_LONG_DOUBLE
//!# Cody-Waite reduction of x such as z=x - n*pi/2 where |z|<=pi/4, result = n mod 8.
#	if MC_TARGET_LONG_DOUBLE_TYPE == MC_TARGET_LONG_DOUBLE_X87
	const long double dp1 = +7.853981554508209228515625000000000000000000000000000000000000000E-01L;
	const long double dp2 = +7.946627356147928367136046290000000000000000000000000000000000000E-09L;
	const long double dp3 = +3.061616997868382943065164830000000000000000000000000000000000000E-17L;
#	else
	const long double dp1 = +7.853981633974483067550664827649598009884357452392578125000000000E-1L;
	const long double dp2 = +2.860594363054915898381331279295066080751126082968574179665700000E-18L;
	const long double dp3 = +2.167952532530945256199261006510837992190580800000000000000000000E-35L;
#	endif

	int64_t r = 0;
	long double w;
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	const int64_t max = INT64_MAX;
#	else
	const int64_t max = mc_cast(const int64_t, 0x7FFFFFFFFFFFFFFF);
#	endif

	*z = 0.0L;
	if (mc_fisvall(x) && mc_fabsl(x) <= 549755813888.0L) {
		if (x == 0.0L) {
			return 0.0L;
		}
		w = mc_floorl(x * MCK_KL(MCK_1_PI_4));
		r = mc_itrunc64l(w - 16.0L * mc_floorl(w / 16.0L));
		if (r < (max - 1)) {
			if (r & 1) {
				++r;
				w = w + 1.0L;
			}
			 r = (mc_cast(uint64_t, r) >> 1) & 7;
			*z = ((x - w * dp1) - w * dp2) - w * dp3;
		} else {
			 r = max;
			*z = MCK_NAN;
		}
	}
	return r;
#	else
	const double y = mc_cast(const double, x);
	double w       = 0.0;
	int64_t r      = mc_rempio2_cw(y, &w);
	*z             = mc_cast(long double, w);
	return r;
#	endif
}

#endif /* !MC_REMPIO2_CW_H */

/* EOF */
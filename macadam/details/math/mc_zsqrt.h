//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zsqrt.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/math/mc_zabs.h>
#include <macadam/details/math/mc_zarg.h>
#include <macadam/details/math/mc_zpolar.h>

#ifndef MC_ZSQRT_H
#define MC_ZSQRT_H

#pragma mark - mc_zsqrt -

MC_TARGET_PROC void mc_zsqrtf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	if (mc_isinf(a_i)) {
		*c_r = MCK_INFP;
		*c_i = a_i;
	} else if (mc_isinf(a_r)) {
		if (mc_isinfp(a_r)) {
			*c_r = a_r;
			*c_i = mc_isnan(a_i) ? a_i : mc_copysignf(0.0f, a_i);
		} else {
			*c_r = mc_isnan(a_i) ? a_i : 0.0f;
			*c_i = mc_copysignf(a_r, a_i);
		}
	} else {
		mc_zpolarf(c_r, c_i
			, mc_sqrtf(mc_zabsf(a_r, a_i))
			, mc_zargf(a_r, a_i) * 0.5f
		);
	}
}

MC_TARGET_PROC void mc_zsqrt(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	if (mc_isinf(a_i)) {
		*c_r = MCK_INFP;
		*c_i = a_i;
	} else if (mc_isinf(a_r)) {
		if (mc_isinfp(a_r)) {
			*c_r = a_r;
			*c_i = mc_isnan(a_i) ? a_i : mc_copysign(0.0, a_i);
		} else {
			*c_r = mc_isnan(a_i) ? a_i : 0.0;
			*c_i = mc_copysign(a_r, a_i);
		}
	} else {
		mc_zpolar(c_r, c_i
			, mc_sqrt(mc_zabs(a_r, a_i))
			, mc_zarg(a_r, a_i) * 0.5
		);
	}
}

MC_TARGET_PROC void mc_zsqrtl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	if (mc_isinf(a_i)) {
		*c_r = MCK_INFP;
		*c_i = a_i;
	} else if (mc_isinf(a_r)) {
		if (mc_isinfp(a_r)) {
			*c_r = a_r;
			*c_i = mc_isnan(a_i) ? a_i : mc_copysignl(0.0L, a_i);
		} else {
			*c_r = mc_isnan(a_i) ? a_i : 0.0L;
			*c_i = mc_copysignl(a_r, a_i);
		}
	} else {
		mc_zpolarl(c_r, c_i
			, mc_sqrtl(mc_zabsl(a_r, a_i))
			, mc_zargl(a_r, a_i) * 0.5L
		);
	}
}

#endif /* !MC_ZSQRT_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zasinh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_zadd.h>
#include <macadam/details/math/mc_zlog.h>
#include <macadam/details/math/mc_zpow.h>
#include <macadam/details/math/mc_zsqrt.h>

#ifndef MC_ZASINH_H
#define MC_ZASINH_H

#pragma mark - mc_zasinh -

MC_TARGET_PROC void mc_zasinhf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	if (mc_isinf(a_r)) {
		if (mc_isnan(a_i)) {
			*c_r = a_r;
			*c_i = a_i;
		} else if (mc_isinf(a_i)) {
			*c_r = a_r;
			*c_i = mc_copysignf(MCK_KF(MCK_PI_4), a_i);
		} else {
			*c_r = a_r;
			*c_i = mc_copysignf(0.0f, a_i);
		}
	} else if (mc_isnan(a_r)) {
		if (mc_isinf(a_i)) {
			*c_r = a_i;
			*c_i = a_r;
		} else if (a_i == 0.0f) {
			*c_r = a_r;
			*c_i = a_i;
		} else {
			*c_r = a_r;
			*c_i = a_r;
		}
	} else if (mc_isinf(a_i)) {
		*c_r = mc_copysignf(a_i, a_r);
		*c_i = mc_copysignf(MCK_KF(MCK_PI_2), a_i);
	} else {
		mc_zpowf(c_r, c_i, a_r, a_i, 2.0f, 0.0f);
		mc_zsqrtf(c_r, c_i, *c_r, *c_i);
		mc_zaddf(c_r, c_i, a_r, a_i, *c_r, *c_i);
		mc_zaddf(c_r, c_i, *c_r, *c_i, 1.0f, 0.0f);
		mc_zlogf(c_r, c_i, *c_r, *c_i);

		*c_r = mc_copysignf(*c_r, a_r);
		*c_i = mc_copysignf(*c_i, a_i);
	}
}

MC_TARGET_PROC void mc_zasinh(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	if (mc_isinf(a_r)) {
		if (mc_isnan(a_i)) {
			*c_r = a_r;
			*c_i = a_i;
		} else if (mc_isinf(a_i)) {
			*c_r = a_r;
			*c_i = mc_copysign(MCK_K(MCK_PI_4), a_i);
		} else {
			*c_r = a_r;
			*c_i = mc_copysign(0.0, a_i);
		}
	} else if (mc_isnan(a_r)) {
		if (mc_isinf(a_i)) {
			*c_r = a_i;
			*c_i = a_r;
		} else if (a_i == 0.0) {
			*c_r = a_r;
			*c_i = a_i;
		} else {
			*c_r = a_r;
			*c_i = a_r;
		}
	} else if (mc_isinf(a_i)) {
		*c_r = mc_copysign(a_i, a_r);
		*c_i = mc_copysign(MCK_K(MCK_PI_2), a_i);
	} else {
		mc_zpow(c_r, c_i, a_r, a_i, 2.0, 0.0);
		mc_zsqrt(c_r, c_i, *c_r, *c_i);
		mc_zadd(c_r, c_i, a_r, a_i, *c_r, *c_i);
		mc_zadd(c_r, c_i, *c_r, *c_i, 1.0, 0.0);
		mc_zlog(c_r, c_i, *c_r, *c_i);

		*c_r = mc_copysign(*c_r, a_r);
		*c_i = mc_copysign(*c_i, a_i);
	}
}

MC_TARGET_PROC void mc_zasinhl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	if (mc_isinf(a_r)) {
		if (mc_isnan(a_i)) {
			*c_r = a_r;
			*c_i = a_i;
		} else if (mc_isinf(a_i)) {
			*c_r = a_r;
			*c_i = mc_copysignl(MCK_KL(MCK_PI_4), a_i);
		} else {
			*c_r = a_r;
			*c_i = mc_copysignl(0.0L, a_i);
		}
	} else if (mc_isnan(a_r)) {
		if (mc_isinf(a_i)) {
			*c_r = a_i;
			*c_i = a_r;
		} else if (a_i == 0.0L) {
			*c_r = a_r;
			*c_i = a_i;
		} else {
			*c_r = a_r;
			*c_i = a_r;
		}
	} else if (mc_isinf(a_i)) {
		*c_r = mc_copysignl(a_i, a_r);
		*c_i = mc_copysignl(MCK_KL(MCK_PI_2), a_i);
	} else {
		mc_zpowl(c_r, c_i, a_r, a_i, 2.0L, 0.0L);
		mc_zsqrtl(c_r, c_i, *c_r, *c_i);
		mc_zaddl(c_r, c_i, a_r, a_i, *c_r, *c_i);
		mc_zaddl(c_r, c_i, *c_r, *c_i, 1.0L, 0.0L);
		mc_zlogl(c_r, c_i, *c_r, *c_i);

		*c_r = mc_copysignl(*c_r, a_r);
		*c_i = mc_copysignl(*c_i, a_i);
	}
}

#endif /* !MC_ZASINH_H */

/* EOF */
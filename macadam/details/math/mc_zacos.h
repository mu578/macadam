//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zacos.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_signbit.h>
#include <macadam/details/math/mc_zadd.h>
#include <macadam/details/math/mc_zlog.h>
#include <macadam/details/math/mc_zpow.h>
#include <macadam/details/math/mc_zsqrt.h>
#include <macadam/details/math/mc_zsub.h>

#ifndef MC_ZACOS_H
#define MC_ZACOS_H

#pragma mark - mc_zacos -

MC_TARGET_PROC void mc_zacosf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	float w;
	if (mc_isinf(a_r)) {
		if (mc_isnan(a_i)) {
			*c_r = a_i;
			*c_i = a_r;
		} else if (mc_isinf(a_i)) {
			if (a_r < 0.0f) {
				*c_r =  MCK_KF(MCK_PI) * 0.75f;
				*c_i = -a_i;
			} else {
				*c_r =  MCK_KF(MCK_PI_4);
				*c_i = -a_i;
			}
		} else if (a_r < 0.0f) {
			*c_r = MCK_KF(MCK_PI);
			*c_i = mc_signbitf(a_i) ? -a_r : a_r;
		} else {
			*c_r = 0.0f;
			*c_i = mc_signbitf(a_i) ? a_r : -a_r;
		}
	} else if (mc_isnan(a_r)) {
		if (mc_isinf(a_i)) {
			*c_r =  a_r;
			*c_i = -a_i;
		} else {
			*c_r = a_r;
			*c_i = a_r;
		}
	} else if (mc_isinf(a_i)) {
		*c_r =  MCK_KF(MCK_PI_2);
		*c_i = -a_i;
	} else if (a_r == 0 && (a_i == 0 || isnan(a_i))) {
		*c_r =  MCK_KF(MCK_PI_2);
		*c_i = -a_i;
	} else {
		mc_zpowf(c_r, c_i, a_r, a_i, 2.0f, 0.0f);
		mc_zsqrtf(c_r, c_i, *c_r, *c_i);
		mc_zaddf(c_r, c_i, a_r, a_i, *c_r, *c_i);
		mc_zsubf(c_r, c_i, *c_r, *c_i, 1.0f, 0.0f);
		mc_zlogf(c_r, c_i, *c_r, *c_i);

		if (mc_signbitf(a_i)) {
			 w   = mc_fabsf(*c_r);
			*c_r = mc_fabsf(*c_i);
			*c_i = w;
		} else {
			 w   =  mc_fabsf(*c_r);
			*c_r =  mc_fabsf(*c_i);
			*c_i = -w;
		}
	}
}

MC_TARGET_PROC void mc_zacos(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	double w;
	if (mc_isinf(a_r)) {
		if (mc_isnan(a_i)) {
			*c_r = a_i;
			*c_i = a_r;
		} else if (mc_isinf(a_i)) {
			if (a_r < 0.0) {
				*c_r =  MCK_K(MCK_PI) * 0.75;
				*c_i = -a_i;
			} else {
				*c_r =  MCK_K(MCK_PI_4);
				*c_i = -a_i;
			}
		} else if (a_r < 0.0) {
			*c_r = MCK_K(MCK_PI);
			*c_i = mc_signbit(a_i) ? -a_r : a_r;
		} else {
			*c_r = 0.0;
			*c_i = mc_signbit(a_i) ? a_r : -a_r;
		}
	} else if (mc_isnan(a_r)) {
		if (mc_isinf(a_i)) {
			*c_r =  a_r;
			*c_i = -a_i;
		} else {
			*c_r = a_r;
			*c_i = a_r;
		}
	} else if (mc_isinf(a_i)) {
		*c_r =  MCK_K(MCK_PI_2);
		*c_i = -a_i;
	} else if (a_r == 0 && (a_i == 0 || isnan(a_i))) {
		*c_r =  MCK_K(MCK_PI_2);
		*c_i = -a_i;
	} else {
		mc_zpow(c_r, c_i, a_r, a_i, 2.0, 0.0);
		mc_zsqrt(c_r, c_i, *c_r, *c_i);
		mc_zadd(c_r, c_i, a_r, a_i, *c_r, *c_i);
		mc_zsub(c_r, c_i, *c_r, *c_i, 1.0, 0.0);
		mc_zlog(c_r, c_i, *c_r, *c_i);

		if (mc_signbit(a_i)) {
			 w   = mc_fabs(*c_r);
			*c_r = mc_fabs(*c_i);
			*c_i = w;
		} else {
			 w   =  mc_fabs(*c_r);
			*c_r =  mc_fabs(*c_i);
			*c_i = -w;
		}
	}
}

MC_TARGET_PROC void mc_zacosl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	long double w;
	if (mc_isinf(a_r)) {
		if (mc_isnan(a_i)) {
			*c_r = a_i;
			*c_i = a_r;
		} else if (mc_isinf(a_i)) {
			if (a_r < 0.0L) {
				*c_r =  MCK_KL(MCK_PI) * 0.75f;
				*c_i = -a_i;
			} else {
				*c_r =  MCK_KL(MCK_PI_4);
				*c_i = -a_i;
			}
		} else if (a_r < 0.0L) {
			*c_r = MCK_KL(MCK_PI);
			*c_i = mc_signbitl(a_i) ? -a_r : a_r;
		} else {
			*c_r = 0.0L;
			*c_i = mc_signbitl(a_i) ? a_r : -a_r;
		}
	} else if (mc_isnan(a_r)) {
		if (mc_isinf(a_i)) {
			*c_r =  a_r;
			*c_i = -a_i;
		} else {
			*c_r = a_r;
			*c_i = a_r;
		}
	} else if (mc_isinf(a_i)) {
		*c_r =  MCK_KL(MCK_PI_2);
		*c_i = -a_i;
	} else if (a_r == 0 && (a_i == 0 || isnan(a_i))) {
		*c_r =  MCK_KL(MCK_PI_2);
		*c_i = -a_i;
	} else {
		mc_zpowl(c_r, c_i, a_r, a_i, 2.0L, 0.0L);
		mc_zsqrtl(c_r, c_i, *c_r, *c_i);
		mc_zaddl(c_r, c_i, a_r, a_i, *c_r, *c_i);
		mc_zsubl(c_r, c_i, *c_r, *c_i, 1.0L, 0.0L);
		mc_zlogl(c_r, c_i, *c_r, *c_i);

		if (mc_signbitl(a_i)) {
			 w   = mc_fabsl(*c_r);
			*c_r = mc_fabsl(*c_i);
			*c_i = w;
		} else {
			 w   =  mc_fabsl(*c_r);
			*c_r =  mc_fabsl(*c_i);
			*c_i = -w;
		}
	}
}

#endif /* !MC_ZACOS_H */

/* EOF */
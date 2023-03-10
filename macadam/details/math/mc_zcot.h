//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zcot.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_cosh.h>
#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_fisnear.h>
#include <macadam/details/math/mc_isfinite.h>
#include <macadam/details/math/mc_sincos.h>
#include <macadam/details/math/mc_sinh.h>
#include <macadam/details/math/mc_zdiv.h>

#ifndef MC_ZCOT_H
#define MC_ZCOT_H

#pragma mark - mc_zcot -

MC_TARGET_PROC void mc_zcotf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = a_r;
	*c_i = a_i;
	if (mc_isfinite(a_r) && mc_isfinite(a_i) && a_r != 0.0f && a_i != 0.0f) {
		const float x = 2.0f * a_r;
		const float y = 2.0f * a_i;
		const float d = mc_coshf(y) - mc_cosf(x);
		if (!mc_fisnearf(d, 0.0f, 2)) {
			const float w = 1.0f / d;
			*c_r          =  mc_sinf (x) * w;
			*c_i          = -mc_sinhf(y) * w;
		}
	}
#	else
	float e, c, s, za_r, za_i, zb_r, zb_i;

	if (mc_isfinite(a_r) && mc_isfinite(a_i) && a_r != 0.0f && a_i != 0.0f) {
		mc_sincosf(2.0f * a_r, &s, &c);

		e    = mc_expf(-2.0f * a_i);
		za_r = -e * s;
		za_i =  e * c;
		zb_r =  za_i - 1.0f;
		za_i =  za_i + 1.0f;
		zb_i = -za_r;

		mc_zdivf(c_r, c_i, za_r, za_i, zb_r, zb_i);
	} else {
		*c_r = a_r;
		*c_i = a_i;
	}
#	endif
}

MC_TARGET_PROC void mc_zcot(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = a_r;
	*c_i = a_i;
	if (mc_isfinite(a_r) && mc_isfinite(a_i) && a_r != 0.0 && a_i != 0.0) {
		const double x = 2.0 * a_r;
		const double y = 2.0 * a_i;
		const double d = mc_cosh(y) - mc_cos(x);
		if (!mc_fisnear(d, 0.0, 2)) {
			const double w = 1.0 / d;
			*c_r           =  mc_sin (x) * w;
			*c_i           = -mc_sinh(y) * w;
		}
	}
#	else
	double e, c, s, za_r, za_i, zb_r, zb_i;

	if (mc_isfinite(a_r) && mc_isfinite(a_i) && a_r != 0.0 && a_i != 0.0) {
		mc_sincos(2.0 * a_r, &s, &c);

		e    = mc_exp(-2.0 * a_i);
		za_r = -e * s;
		za_i =  e * c;
		zb_r =  za_i - 1.0;
		za_i =  za_i + 1.0;
		zb_i = -za_r;

		mc_zdiv(c_r, c_i, za_r, za_i, zb_r, zb_i);
	} else {
		*c_r = a_r;
		*c_i = a_i;
	}
#	endif
}

MC_TARGET_PROC void mc_zcotl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = a_r;
	*c_i = a_i;
	if (mc_isfinite(a_r) && mc_isfinite(a_i) && a_r != 0.0L && a_i != 0.0L) {
		const long double x = 2.0L * a_r;
		const long double y = 2.0L * a_i;
		const long double d = mc_coshl(y) - mc_cosl(x);
		if (!mc_fisnearl(d, 0.0L, 2)) {
			const long double w = 1.0L / d;
			*c_r                =  mc_sinl (x) * w;
			*c_i                = -mc_sinhl(y) * w;
		}
	}
#	else
	long double e, c, s, za_r, za_i, zb_r, zb_i;

	if (mc_isfinite(a_r) && mc_isfinite(a_i) && a_r != 0.0L && a_i != 0.0L) {
		mc_sincosl(2.0L * a_r, &s, &c);

		e    = mc_expl(-2.0L * a_i);
		za_r = -e * s;
		za_i =  e * c;
		zb_r =  za_i - 1.0L;
		za_i =  za_i + 1.0L;
		zb_i = -za_r;

		mc_zdivl(c_r, c_i, za_r, za_i, zb_r, zb_i);
	} else {
		*c_r = a_r;
		*c_i = a_i;
	}
#	endif
}

#endif /* !MC_ZCOT_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zexpm1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_expm1.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sin.h>
#include <macadam/details/math/mc_zexp.h>
#include <macadam/details/math/mc_zsub.h>

#ifndef MC_ZEXPM1_H
#define MC_ZEXPM1_H

#pragma mark - mc_zexpm1 -

MC_TARGET_PROC void mc_zexpm1f(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	if ((mc_fabsf(a_r) >= 1.0f) || (mc_fabsf(a_i) >= 1.0f)) {
		mc_zexpf(c_r, c_i, a_r, a_i);
		mc_zsubf(c_r, c_i, *c_r, *c_i, 1.0f, 0.0f);
	} else {
		const float zm1 = mc_expm1f(a_r);
		const float zp1 = 1.0f + zm1;
		*c_r            = zm1 - 2.0f * zp1 * mc_raise2f(mc_sinf(0.5f * a_i));
		*c_i            = zp1 * mc_sinf(a_i);
	}
}

MC_TARGET_PROC void mc_zexpm1(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	if ((mc_fabs(a_r) >= 1.0) || (mc_fabs(a_i) >= 1.0)) {
		mc_zexp(c_r, c_i, a_r, a_i);
		mc_zsub(c_r, c_i, *c_r, *c_i, 1.0, 0.0);
	} else {
		const double zm1 = mc_expm1(a_r);
		const double zp1 = 1.0 + zm1;
		*c_r             = zm1 - 2.0 * zp1 * mc_raise2(mc_sin(0.5 * a_i));
		*c_i             = zp1 * mc_sin(a_i);
	}
}

MC_TARGET_PROC void mc_zexpm1l(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	if ((mc_fabsl(a_r) >= 1.0L) || (mc_fabsl(a_i) >= 1.0L)) {
		mc_zexpl(c_r, c_i, a_r, a_i);
		mc_zsubl(c_r, c_i, *c_r, *c_i, 1.0L, 0.0L);
	} else {
		const long double zm1 = mc_expm1l(a_r);
		const long double zp1 = 1.0L + zm1;
		*c_r                  = zm1 - 2.0L * zp1 * mc_raise2l(mc_sinl(0.5L * a_i));
		*c_i                  = zp1 * mc_sinl(a_i);
	}
}

#endif /* !MC_ZEXPM1_H */

/* EOF */
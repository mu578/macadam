//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zlog1p.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_atan2.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_log1p.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_zadd.h>
#include <macadam/details/math/mc_zlog.h>

#ifndef MC_ZLOG1P_H
#define MC_ZLOG1P_H

#pragma mark - mc_zlog1p -

MC_TARGET_PROC void mc_zlog1pf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	if ((mc_fabsf(a_r) >= 1.0f) || (mc_fabsf(a_i) >= 1.0f)) {
		mc_zaddf(c_r, c_i, 1.0f, 0.0f, a_r, a_i);
		mc_zlogf(c_r, c_i, *c_r, *c_i);
	} else {
		const float zp1 = 1.0f + a_r;
		*c_r            = mc_log1pf(a_r) + 0.5f * mc_log1pf(mc_raise2f(a_i / zp1));
		*c_i            = mc_atan2f(a_i, zp1);
	}
}

MC_TARGET_PROC void mc_zlog1p(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	if ((mc_fabs(a_r) >= 1.0) || (mc_fabs(a_i) >= 1.0)) {
		mc_zadd(c_r, c_i, 1.0, 0.0, a_r, a_i);
		mc_zlog(c_r, c_i, *c_r, *c_i);
	} else {
		const double zp1 = 1.0 + a_r;
		*c_r             = mc_log1p(a_r) + 0.5 * mc_log1p(mc_raise2(a_i / zp1));
		*c_i             = mc_atan2(a_i, zp1);
	}
}

MC_TARGET_PROC void mc_zlog1pl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	if ((mc_fabsl(a_r) >= 1.0L) || (mc_fabsl(a_i) >= 1.0L)) {
		mc_zaddl(c_r, c_i, 1.0L, 0.0L, a_r, a_i);
		mc_zlogl(c_r, c_i, *c_r, *c_i);
	} else {
		const long double zp1 = 1.0L + a_r;
		*c_r                  = mc_log1pl(a_r) + 0.5L * mc_log1pl(mc_raise2l(a_i / zp1));
		*c_i                  = mc_atan2l(a_i, zp1);
	}
}

#endif /* !MC_ZLOG1P_H */

/* EOF */
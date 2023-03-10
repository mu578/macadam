//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zpochhammer.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zmul.h>

#ifndef MC_ZPOCHHAMMER_H
#define MC_ZPOCHHAMMER_H

#pragma mark - mc_zpochhammer -

MC_TARGET_PROC void zpochhammerf(
	  float *     r_r
	, float *     r_i
	, const float x_r
	, const float x_i
	, const int   n
) {
	int i = 1;
	*r_r  = 1.0f;
	*r_i  = 0.0f;
	if (n > 0 && n < MCLIMITS_IMAX) {
		for (; i <= n; i++) {
			mc_zmulf(r_r, r_i, *r_r, *r_i, (x_r + mc_cast(float, i)) - 1.0f, x_i);
		}
	}
}

MC_TARGET_PROC void mc_zpochhammer(
	  double *     r_r
	, double *     r_i
	, const double x_r
	, const double x_i
	, const int    n
) {
	int i = 1;
	*r_r  = 1.0;
	*r_i  = 0.0;
	if (n > 0 && n < MCLIMITS_IMAX) {
		for (; i <= n; i++) {
			mc_zmul(r_r, r_i, *r_r, *r_i, (x_r + mc_cast(double, i)) - 1.0, x_i);
		}
	}
}

MC_TARGET_PROC void mc_zpochhammerl(
	  long double *     r_r
	, long double *     r_i
	, const long double x_r
	, const long double x_i
	, const int         n
) {
	int i = 1;
	*r_r  = 1.0L;
	*r_i  = 0.0L;
	if (n > 0 && n < MCLIMITS_IMAX) {
		for (; i <= n; i++) {
			mc_zmull(r_r, r_i, *r_r, *r_i, (x_r + mc_cast(long double, i)) - 1.0L, x_i);
		}
	}
}

#endif /* !MC_ZPOCHHAMMER_H */

/* EOF */
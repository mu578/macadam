//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zrecip.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zabs.h>

#ifndef MC_ZRECIP_H
#define MC_ZRECIP_H

#pragma mark - mc_zrecip -

MC_TARGET_PROC void mc_zrecipf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
//!# Similar to @zinv for internal use.
	const float d = mc_zabsf(a_r, a_i);
	const float r = d != 0.0f ? 1.0f / d : MCK_NAN;
	if (!mc_isnan(r)) {
		*c_r = +(a_r * r) * r;
		*c_i = -(a_i * r) * r;
	} else {
		*c_r = MCK_NAN;
		*c_i = MCK_NAN;
	}
}

MC_TARGET_PROC void mc_zrecip(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
//!# Similar to @zinv for internal use.
	const double d = mc_zabs(a_r, a_i);
	const double r = d != 0.0 ? 1.0 / d : MCK_NAN;
	if (!mc_isnan(r)) {
		*c_r = +(a_r * r) * r;
		*c_i = -(a_i * r) * r;
	} else {
		*c_r = MCK_NAN;
		*c_i = MCK_NAN;
	}
}

MC_TARGET_PROC void mc_zrecipl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
//!# Similar to @zinv for internal use.
	const long double d = mc_zabsl(a_r, a_i);
	const long double r = d != 0.0L ? 1.0L / d : MCK_NAN;
	if (!mc_isnan(r)) {
		*c_r = +(a_r * r) * r;
		*c_i = -(a_i * r) * r;
	} else {
		*c_r = MCK_NAN;
		*c_i = MCK_NAN;
	}
}

#endif /* !MC_ZRECIP_H */

/* EOF */
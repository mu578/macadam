//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zinv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_zdiv.h>

#ifndef MC_ZINV_H
#define MC_ZINV_H

#pragma mark - mc_zinv -

MC_TARGET_PROC void mc_zinvf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	const float d = mc_zabsf(a_r, a_i);
	const float r = d != 0.0f ? 1.0f / d : MCK_NAN;
	if (!mc_isnan(r)) {
		*c_r = +(a_r * r) * r;
		*c_i = -(a_i * r) * r;
	} else {
		*c_r = 0.0f;
		*c_i = -mc_copysignf(0.0f, a_i);
	}
}

MC_TARGET_PROC void mc_zinv(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	const double d = mc_zabs(a_r, a_i);
	const double r = d != 0.0 ? 1.0 / d : MCK_NAN;
	if (!mc_isnan(r)) {
		*c_r = +(a_r * r) * r;
		*c_i = -(a_i * r) * r;
	} else {
		*c_r = 0.0;
		*c_i = -mc_copysign(0.0, a_i);
	}
}

MC_TARGET_PROC void mc_zinvl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	const long double d = mc_zabsl(a_r, a_i);
	const long double r = d != 0.0L ? 1.0L / d : MCK_NAN;
	if (!mc_isnan(r)) {
		*c_r = +(a_r * r) * r;
		*c_i = -(a_i * r) * r;
	} else {
		*c_r = 0.0L;
		*c_i = -mc_copysignl(0.0L, a_i);
	}
}

#endif /* !MC_ZINV_H */

/* EOF */
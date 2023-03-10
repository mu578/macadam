//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zexp2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zpow.h>

#ifndef MC_ZEXP2_H
#define MC_ZEXP2_H

#pragma mark - mc_zexp2 -

MC_TARGET_PROC void mc_zexp2f(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	*c_r = 0.0f;
	*c_i = 0.0f;
	mc_zpowf(c_r, c_i, 2.0f, 0.0f, a_r, a_i);
}

MC_TARGET_PROC void mc_zexp2(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	mc_zpow(c_r, c_i, 2.0, 0.0, a_r, a_i);
}

MC_TARGET_PROC void mc_zexp2l(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	mc_zpowl(c_r, c_i, 2.0L, 0.0L, a_r, a_i);
}

#endif /* !MC_ZEXP2_H */

/* EOF */
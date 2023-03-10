//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zfpow.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zpow.h>

#ifndef MC_ZFPOW_H
#define MC_ZFPOW_H

#pragma mark - mc_zfpow -

MC_TARGET_PROC void mc_zfpowf(float * c_r, float * c_i
	, const float a_r, const float a_i
	, const float b
) {
	mc_zpowf(c_r, c_i, a_r, a_i, b, 0.0f);
}

MC_TARGET_PROC void mc_zfpow(double * c_r, double * c_i
	, const double a_r, const double a_i
	, const double b
) {
	mc_zpow(c_r, c_i, a_r, a_i, b, 0.0);
}

MC_TARGET_PROC void mc_zfpowl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
	, const long double b
) {
	mc_zpowl(c_r, c_i, a_r, a_i, b, 0.0L);
}

#endif /* !MC_ZFPOW_H */

/* EOF */
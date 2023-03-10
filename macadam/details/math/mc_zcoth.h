//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zcoth.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zrecip.h>
#include <macadam/details/math/mc_ztanh.h>

#ifndef MC_ZCOTH_H
#define MC_ZCOTH_H

#pragma mark - mc_zcoth -

MC_TARGET_PROC void mc_zcothf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	mc_ztanhf(c_r, c_i, a_r, a_i);
	mc_zrecipf(c_r, c_i, *c_r, *c_i);
}

MC_TARGET_PROC void mc_zcoth(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	mc_ztanh(c_r, c_i, a_r, a_i);
	mc_zrecip(c_r, c_i, *c_r, *c_i);
}

MC_TARGET_PROC void mc_zcothl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	mc_ztanhl(c_r, c_i, a_r, a_i);
	mc_zrecipl(c_r, c_i, *c_r, *c_i);
}

#endif /* !MC_ZCOTH_H */

/* EOF */
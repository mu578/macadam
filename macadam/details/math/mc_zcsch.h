//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zcsch.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zsinh.h>

#ifndef MC_ZCSCH_H
#define MC_ZCSCH_H

#pragma mark - mc_zcsch -

MC_TARGET_PROC void mc_zcschf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	mc_zsinhf(c_r, c_i, a_r, a_i);
	mc_zrecipf(c_r, c_i, *c_r, *c_i);
}

MC_TARGET_PROC void mc_zcsch(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	mc_zsinh(c_r, c_i, a_r, a_i);
	mc_zrecip(c_r, c_i, *c_r, *c_i);
}

MC_TARGET_PROC void mc_zcschl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	mc_zsinhl(c_r, c_i, a_r, a_i);
	mc_zrecipl(c_r, c_i, *c_r, *c_i);
}

#endif /* !MC_ZCSCH_H */

/* EOF */
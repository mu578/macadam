//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zcsc.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zsin.h>

#ifndef MC_ZCSC_H
#define MC_ZCSC_H

#pragma mark - mc_zcsc -

MC_TARGET_PROC void mc_zcscf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	mc_zsinf(c_r, c_i, a_r, a_i);
	mc_zrecipf(c_r, c_i, *c_r, *c_i);
}

MC_TARGET_PROC void mc_zcsc(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	mc_zsin(c_r, c_i, a_r, a_i);
	mc_zrecip(c_r, c_i, *c_r, *c_i);
}

MC_TARGET_PROC void mc_zcscl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	mc_zsinl(c_r, c_i, a_r, a_i);
	mc_zrecipl(c_r, c_i, *c_r, *c_i);
}

#endif /* !MC_ZCSC_H */

/* EOF */
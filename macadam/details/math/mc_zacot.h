//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zacot.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zatan.h>
#include <macadam/details/math/mc_zrecip.h>

#ifndef MC_ZACOT_H
#define MC_ZACOT_H

#pragma mark - mc_zacot -

MC_TARGET_PROC void mc_zacotf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	mc_zrecipf(c_r, c_i, a_r, a_i);
	mc_zatanf(c_r, c_i, *c_r, *c_i);
}

MC_TARGET_PROC void mc_zacot(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	mc_zrecip(c_r, c_i, a_r, a_i);
	mc_zatan(c_r, c_i, *c_r, *c_i);
}

MC_TARGET_PROC void mc_zacotl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	mc_zrecipl(c_r, c_i, a_r, a_i);
	mc_zatanl(c_r, c_i, *c_r, *c_i);
}

#endif /* !MC_ZACOT_H */

/* EOF */
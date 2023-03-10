//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zacoth.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zatanh.h>
#include <macadam/details/math/mc_zrecip.h>

#ifndef MC_ZACOTH_H
#define MC_ZACOTH_H

#pragma mark - mc_zacoth -

MC_TARGET_PROC void mc_zacothf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	mc_zrecipf(c_r, c_i, a_r, a_i);
	mc_zatanhf(c_r, c_i, *c_r, *c_i);
}

MC_TARGET_PROC void mc_zacoth(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	mc_zrecip(c_r, c_i, a_r, a_i);
	mc_zatanh(c_r, c_i, *c_r, *c_i);
}

MC_TARGET_PROC void mc_zacothl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	mc_zrecipl(c_r, c_i, a_r, a_i);
	mc_zatanhl(c_r, c_i, *c_r, *c_i);
}

#endif /* !MC_ZACOTH_H */

/* EOF */
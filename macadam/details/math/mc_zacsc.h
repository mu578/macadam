//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zacsc.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zasin.h>
#include <macadam/details/math/mc_zrecip.h>

#ifndef MC_ZACSC_H
#define MC_ZACSC_H

#pragma mark - mc_zacsc -

MC_TARGET_PROC void mc_zacscf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	mc_zrecipf(c_r, c_i, a_r, a_i);
	mc_zasinf(c_r, c_i, *c_r, *c_i);
}

MC_TARGET_PROC void mc_zacsc(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	mc_zrecip(c_r, c_i, a_r, a_i);
	mc_zasin(c_r, c_i, *c_r, *c_i);
}

MC_TARGET_PROC void mc_zacscl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	mc_zrecipl(c_r, c_i, a_r, a_i);
	mc_zasinl(c_r, c_i, *c_r, *c_i);
}

#endif /* !MC_ZACSC_H */

/* EOF */
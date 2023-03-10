//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zatan.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zatanh.h>

#ifndef MC_ZATAN_H
#define MC_ZATAN_H

#pragma mark - mc_zatan -

MC_TARGET_PROC void mc_zatanf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	float w;
	mc_zatanhf(c_r, c_i, -a_i, a_r);
	 w   = *c_r;
	*c_r = *c_i;
	*c_i = -w;
}

MC_TARGET_PROC void mc_zatan(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	double w;
	mc_zatanh(c_r, c_i, -a_i, a_r);
	 w   = *c_r;
	*c_r = *c_i;
	*c_i = -w;
}

MC_TARGET_PROC void mc_zatanl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	long double w;
	mc_zatanhl(c_r, c_i, -a_i, a_r);
	 w   = *c_r;
	*c_r = *c_i;
	*c_i = -w;
}

#endif /* !MC_ZATAN_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ztan.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_ztanh.h>

#ifndef MC_ZTAN_H
#define MC_ZTAN_H

#pragma mark - mc_ztan -

MC_TARGET_PROC void mc_ztanf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	float w;
	mc_ztanhf(c_r, c_i, -a_i, a_r);
	 w   = *c_r;
	*c_r = *c_i;
	*c_i = -w;
}

MC_TARGET_PROC void mc_ztan(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	double w;
	mc_ztanh(c_r, c_i, -a_i, a_r);
	 w   = *c_r;
	*c_r = *c_i;
	*c_i = -w;
}

MC_TARGET_PROC void mc_ztanl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	long double w;
	mc_ztanhl(c_r, c_i, -a_i, a_r);
	 w   = *c_r;
	*c_r = *c_i;
	*c_i = -w;
}

#endif /* !MC_ZTAN_H */

/* EOF */
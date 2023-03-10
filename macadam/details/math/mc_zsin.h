//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zsin.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_cos.h>
#include <macadam/details/math/mc_cosh.h>
#include <macadam/details/math/mc_sin.h>
#include <macadam/details/math/mc_sinh.h>
#include <macadam/details/math/mc_zsinh.h>

#ifndef MC_ZSIN_H
#define MC_ZSIN_H

#pragma mark - mc_zsin -

MC_TARGET_PROC void mc_zsinf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = mc_sinf(a_r) * mc_coshf(a_i);
	*c_i = mc_cosf(a_r) * mc_sinhf(a_i);
#	else
	float w;
	mc_zsinhf(c_r, c_i, -a_i, a_r);
	 w   = *c_r;
	*c_r = *c_i;
	*c_i = -w;
#	endif
}

MC_TARGET_PROC void mc_zsin(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = mc_sin(a_r) * mc_cosh(a_i);
	*c_i = mc_cos(a_r) * mc_sinh(a_i);
#	else
	double w;
	mc_zsinh(c_r, c_i, -a_i, a_r);
	 w   = *c_r;
	*c_r = *c_i;
	*c_i = -w;
#	endif
}

MC_TARGET_PROC void mc_zsinl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = mc_sinl(a_r) * mc_coshl(a_i);
	*c_i = mc_cosl(a_r) * mc_sinhl(a_i);
#	else
	long double w;
	mc_zsinhl(c_r, c_i, -a_i, a_r);
	 w   = *c_r;
	*c_r = *c_i;
	*c_i = -w;
#	endif
}

#endif /* !MC_ZSIN_H */

/* EOF */
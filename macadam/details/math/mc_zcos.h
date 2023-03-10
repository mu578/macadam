//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zcos.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_cos.h>
#include <macadam/details/math/mc_cosh.h>
#include <macadam/details/math/mc_sin.h>
#include <macadam/details/math/mc_sinh.h>
#include <macadam/details/math/mc_zcosh.h>

#ifndef MC_ZCOS_H
#define MC_ZCOS_H

#pragma mark - mc_zcos -

MC_TARGET_PROC void mc_zcosf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r =  mc_cosf(a_r) * mc_coshf(a_i);
	*c_i = -mc_sinf(a_r) * mc_sinhf(a_i);
#	else
	mc_zcoshf(c_r, c_i, -a_i, a_r);
#	endif
}

MC_TARGET_PROC void mc_zcos(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r =  mc_cos(a_r) * mc_cosh(a_i);
	*c_i = -mc_sin(a_r) * mc_sinh(a_i);
#	else
	mc_zcosh(c_r, c_i, -a_i, a_r);
#	endif
}

MC_TARGET_PROC void mc_zcosl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r =  mc_cosl(a_r) * mc_coshl(a_i);
	*c_i = -mc_sinl(a_r) * mc_sinhl(a_i);
#	else
	mc_zcoshl(c_r, c_i, -a_i, a_r);
#	endif
}

#endif /* !MC_ZCOS_H */

/* EOF */
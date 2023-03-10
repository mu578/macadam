//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zsinh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_cos.h>
#include <macadam/details/math/mc_cosh.h>
#include <macadam/details/math/mc_isfinite.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_sin.h>
#include <macadam/details/math/mc_sinh.h>

#ifndef MC_ZSINH_H
#define MC_ZSINH_H

#pragma mark - mc_zsinh -

MC_TARGET_PROC void mc_zsinhf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = mc_sinhf(a_r) * mc_cosf(a_i);
	*c_i = mc_coshf(a_r) * mc_sinf(a_i);
#	else
	if ((mc_isinf(a_r) || a_r == 0.0f) && !mc_isfinite(a_i)) {
		*c_r = a_r;
		*c_i = MCK_NAN;
	} else if (a_i == 0.0f && !mc_isfinite(a_r)) {
		*c_r = a_r;
		*c_i = a_i;
	} else {
		*c_r = mc_sinhf(a_r) * mc_cosf(a_i);
		*c_i = mc_coshf(a_r) * mc_sinf(a_i);
	}
#	endif
}

MC_TARGET_PROC void mc_zsinh(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = mc_sinh(a_r) * mc_cos(a_i);
	*c_i = mc_cosh(a_r) * mc_sin(a_i);
#	else
	if ((mc_isinf(a_r) || a_r == 0.0) && !mc_isfinite(a_i)) {
		*c_r = a_r;
		*c_i = MCK_NAN;
	} else if (a_i == 0.0 && !mc_isfinite(a_r)) {
		*c_r = a_r;
		*c_i = a_i;
	} else {
		*c_r = mc_sinh(a_r) * mc_cos(a_i);
		*c_i = mc_cosh(a_r) * mc_sin(a_i);
	}
#	endif
}

MC_TARGET_PROC void mc_zsinhl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = mc_sinhl(a_r) * mc_cosl(a_i);
	*c_i = mc_coshl(a_r) * mc_sinl(a_i);
#	else
	if ((mc_isinf(a_r) || a_r == 0.0L) && !mc_isfinite(a_i)) {
		*c_r = a_r;
		*c_i = MCK_NAN;
	} else if (a_i == 0.0L && !mc_isfinite(a_r)) {
		*c_r = a_r;
		*c_i = a_i;
	} else {
		*c_r = mc_sinhl(a_r) * mc_cosl(a_i);
		*c_i = mc_coshl(a_r) * mc_sinl(a_i);
	}
#	endif
}

#endif /* !MC_ZSINH_H */

/* EOF */
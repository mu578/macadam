//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zcosh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_cos.h>
#include <macadam/details/math/mc_cosh.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_isfinite.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_sin.h>
#include <macadam/details/math/mc_sinh.h>

#ifndef MC_ZCOSH_H
#define MC_ZCOSH_H

#pragma mark - mc_zcosh -

MC_TARGET_PROC void mc_zcoshf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = mc_coshf(a_r) * mc_cosf(a_i);
	*c_i = mc_sinhf(a_r) * mc_sinf(a_i);
#	else
	if (mc_isinf(a_r) && !mc_isfinite(a_i)) {
		*c_r = mc_fabsf(a_r);
		*c_i = MCK_NAN;
	} else if (a_r == 0.0f && !mc_isfinite(a_i)) {
		*c_r = MCK_NAN;
		*c_i = a_r;
	} else if (a_r == 0.0f && a_i == 0.0f) {
		*c_r = 1.0f;
		*c_i = a_i;
	} else if (a_i == 0.0f && !mc_isfinite(a_r)) {
		*c_r = mc_fabsf(a_r);
		*c_i = a_i;
	} else {
		*c_r = mc_coshf(a_r) * mc_cosf(a_i);
		*c_i = mc_sinhf(a_r) * mc_sinf(a_i);
	}
#	endif
}

MC_TARGET_PROC void mc_zcosh(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = mc_cosh(a_r) * mc_cos(a_i);
	*c_i = mc_sinh(a_r) * mc_sin(a_i);
#	else
	if (mc_isinf(a_r) && !mc_isfinite(a_i)) {
		*c_r = mc_fabs(a_r);
		*c_i = MCK_NAN;
	} else if (a_r == 0.0 && !mc_isfinite(a_i)) {
		*c_r = MCK_NAN;
		*c_i = a_r;
	} else if (a_r == 0.0 && a_i == 0.0) {
		*c_r = 1.0;
		*c_i = a_i;
	} else if (a_i == 0.0 && !mc_isfinite(a_r)) {
		*c_r = mc_fabs(a_r);
		*c_i = a_i;
	} else {
		*c_r = mc_cosh(a_r) * mc_cos(a_i);
		*c_i = mc_sinh(a_r) * mc_sin(a_i);
	}
#	endif
}

MC_TARGET_PROC void mc_zcoshl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
#	if MC_TARGET_EMBEDDED
	*c_r = mc_coshl(a_r) * mc_cosl(a_i);
	*c_i = mc_sinhl(a_r) * mc_sinl(a_i);
#	else
	if (mc_isinf(a_r) && !mc_isfinite(a_i)) {
		*c_r = mc_fabsl(a_r);
		*c_i = MCK_NAN;
	} else if (a_r == 0.0L && !mc_isfinite(a_i)) {
		*c_r = MCK_NAN;
		*c_i = a_r;
	} else if (a_r == 0.0L && a_i == 0.0L) {
		*c_r = 1.0L;
		*c_i = a_i;
	} else if (a_i == 0.0L && !mc_isfinite(a_r)) {
		*c_r = mc_fabsl(a_r);
		*c_i = a_i;
	} else {
		*c_r = mc_coshl(a_r) * mc_cosl(a_i);
		*c_i = mc_sinhl(a_r) * mc_sinl(a_i);
	}
#	endif
}

#endif /* !MC_ZCOSH_H */

/* EOF */
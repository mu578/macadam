//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zexp.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_isfinite.h>

#ifndef MC_ZEXP_H
#define MC_ZEXP_H

#pragma mark - mc_zexp -

MC_TARGET_PROC void mc_zexpf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
#	if !MC_TARGET_EMBEDDED
	*c_i = a_i;
	if (mc_isinf(a_r)) {
		if (a_r < 0.0f) {
			if (!mc_isfinite(*c_i)) {
				*c_i = 1.0f;
			}
		} else if (*c_i == 0.0f || !mc_isfinite(*c_i)) {
			if (mc_isinf(*c_i)) {
				*c_i = MCK_NAN;
			}
			*c_r = a_r;
			return;
		}
	} else if (mc_isnan(a_r) && a_i == 0.0f) {
		*c_r = a_r;
		*c_i = a_i;
		return;
	}
#	else
	if (a_r == 0.0f && a_i == 0.0f) {
		*c_r = 1.0f;
		*c_i = 0.0f;
		return;
	}
#	endif
	const float x = mc_expf(a_r);
	*c_r          = x * mc_cosf(a_i);
	*c_i          = x * mc_sinf(a_i);
}

MC_TARGET_PROC void mc_zexp(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
#	if !MC_TARGET_EMBEDDED
	*c_i = a_i;
	if (mc_isinf(a_r)) {
		if (a_r < 0.0) {
			if (!mc_isfinite(*c_i)) {
				*c_i = 1.0;
			}
		} else if (*c_i == 0.0 || !mc_isfinite(*c_i)) {
			if (mc_isinf(*c_i)) {
				*c_i = MCK_NAN;
			}
			*c_r = a_r;
			return;
		}
	} else if (mc_isnan(a_r) && a_i == 0.0) {
		*c_r = a_r;
		*c_i = a_i;
		return;
	}
#	else
	if (a_r == 0.0 && a_i == 0.0) {
		*c_r = 1.0;
		*c_i = 0.0;
		return;
	}
#	endif
	const double x = mc_exp(a_r);
	*c_r           = x * mc_cos(a_i);
	*c_i           = x * mc_sin(a_i);
}

MC_TARGET_PROC void mc_zexpl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
#	if !MC_TARGET_EMBEDDED
	*c_i = a_i;
	if (mc_isinf(a_r)) {
		if (a_r < 0.0L) {
			if (!mc_isfinite(*c_i)) {
				*c_i = 1.0L;
			}
		} else if (*c_i == 0.0L || !mc_isfinite(*c_i)) {
			if (mc_isinf(*c_i)) {
				*c_i = MCK_NAN;
			}
			*c_r = a_r;
			return;
		}
	} else if (mc_isnan(a_r) && a_i == 0.0L) {
		*c_r = a_r;
		*c_i = a_i;
		return;
	}
#	else
	if (a_r == 0.0L && a_i == 0.0L) {
		*c_r = 1.0L;
		*c_i = 0.0L;
		return;
	}
#	endif
	const long double x = mc_expl(a_r);
	*c_r                = x * mc_cosl(a_i);
	*c_i                = x * mc_sinl(a_i);
}

#endif /* !MC_ZEXP_H */

/* EOF */
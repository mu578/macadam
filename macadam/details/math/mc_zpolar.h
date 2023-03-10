//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zpolar.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_cos.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_signbit.h>
#include <macadam/details/math/mc_sin.h>

#ifndef MC_ZPOLAR_H
#define MC_ZPOLAR_H

#pragma mark - mc_zpolar -

MC_TARGET_PROC void mc_zpolarf(float * z_r, float * z_i
	, float rho
	, float theta
) {
	if (mc_isnan(rho) || mc_signbitf(rho) != 0) {
		*z_r = MCK_NAN;
		*z_i = MCK_NAN;
	} else if (mc_isnan(theta)) {
		if (mc_isinf(rho)) {
			*z_r = rho;
			*z_i = theta;
		} else {
			*z_r = theta;
			*z_i = theta;
		}
	} else if (mc_isinf(theta)) {
		if (mc_isinf(rho)) {
			*z_r = rho;
			*z_i = MCK_NAN;
		} else {
			*z_r = MCK_NAN;
			*z_i = theta;
		}
	} else {
		*z_r = rho * mc_cosf(theta);
		if (mc_isnan(*z_r)) {
			*z_r = 0.0f;
		}
		*z_i = rho * mc_sinf(theta);
		if (mc_isnan(*z_i)) {
			*z_i = 0.0f;
		}
	}
}

MC_TARGET_PROC void mc_zpolar(double * z_r, double * z_i
	, double rho
	, double theta
) {
	if (mc_isnan(rho) || mc_signbit(rho) != 0) {
		*z_r = MCK_NAN;
		*z_i = MCK_NAN;
	} else if (mc_isnan(theta)) {
		if (mc_isinf(rho)) {
			*z_r = rho;
			*z_i = theta;
		} else {
			*z_r = theta;
			*z_i = theta;
		}
	} else if (mc_isinf(theta)) {
		if (mc_isinf(rho)) {
			*z_r = rho;
			*z_i = MCK_NAN;
		} else {
			*z_r = MCK_NAN;
			*z_i = theta;
		}
	} else {
		*z_r = rho * mc_cos(theta);
		if (mc_isnan(*z_r)) {
			*z_r = 0.0;
		}
		*z_i = rho * mc_sin(theta);
		if (mc_isnan(*z_i)) {
			*z_i = 0.0;
		}
	}
}

MC_TARGET_PROC void mc_zpolarl(long double * z_r, long double * z_i
	, long double rho
	, long double theta
) {
	if (mc_isnan(rho) || mc_signbitl(rho) != 0) {
		*z_r = MCK_NAN;
		*z_i = MCK_NAN;
	} else if (mc_isnan(theta)) {
		if (mc_isinf(rho)) {
			*z_r = rho;
			*z_i = theta;
		} else {
			*z_r = theta;
			*z_i = theta;
		}
	} else if (mc_isinf(theta)) {
		if (mc_isinf(rho)) {
			*z_r = rho;
			*z_i = MCK_NAN;
		} else {
			*z_r = MCK_NAN;
			*z_i = theta;
		}
	} else {
		*z_r = rho * mc_cosl(theta);
		if (mc_isnan(*z_r)) {
			*z_r = 0.0L;
		}
		*z_i = rho * mc_sinl(theta);
		if (mc_isnan(*z_i)) {
			*z_i = 0.0L;
		}
	}
}

#endif /* !MC_ZPOLAR_H */

/* EOF */
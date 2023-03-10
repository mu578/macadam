//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ztanh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_cos.h>
#include <macadam/details/math/mc_cosh.h>
#include <macadam/details/math/mc_isfinite.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_sin.h>
#include <macadam/details/math/mc_sinh.h>

#ifndef MC_ZTANH_H
#define MC_ZTANH_H

#pragma mark - mc_ztanh -

MC_TARGET_PROC void mc_ztanhf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
#	if MC_TARGET_EMBEDDED
	const float z2r = 2.0f * a_r;
	const float z2i = 2.0f * a_i;
	const float zd  = mc_coshf(z2r) + mc_cosf(z2i);
	*c_r            = mc_sinhf(z2r) / zd;
	*c_i            = mc_sinf(z2i) / zd;
#	else
	if (mc_isinf(a_r)) {
		if (!mc_isfinite(a_i)) {
			*c_r = 1.0f;
			*c_i = 0.0f;
		} else {
			*c_r = 1.0f;
			*c_i = mc_copysignf(0.0f, mc_sinf(2.0f * a_i));
		}
	} else if (mc_isnan(a_r) && a_i == 0.0f) {
		*c_r = a_r;
		*c_i = a_i;
	} else {
		const float z2r   = 2.0f * a_r;
		const float z2i   = 2.0f * a_i;
		const float zd    = mc_coshf(z2r) + mc_cosf(z2i);
		const float z2rsh = mc_sinhf(z2r);
		if (mc_isinf(z2rsh) && mc_isinf(zd)) {
			*c_r = z2rsh > 0.0f ? 1.0f : -1.0f;
			*c_i = z2i   > 0.0f ? 0.0f : -0.0f;
		} else {
			*c_r = z2rsh / zd;
			*c_i = mc_sinf(z2i) / zd;
		}
	}
#	endif
}

MC_TARGET_PROC void mc_ztanh(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
#	if MC_TARGET_EMBEDDED
	const double z2r = 2.0 * a_r;
	const double z2i = 2.0 * a_i;
	const double zd  = mc_cosh(z2r) + mc_cos(z2i);
	*c_r             = mc_sinh(z2r) / zd;
	*c_i             = mc_sin(z2i) / zd;
#	else
	if (mc_isinf(a_r)) {
		if (!mc_isfinite(a_i)) {
			*c_r = 1.0;
			*c_i = 0.0;
		} else {
			*c_r = 1.0;
			*c_i = mc_copysign(0.0, mc_sin(2.0 * a_i));
		}
	} else if (mc_isnan(a_r) && a_i == 0.0) {
		*c_r = a_r;
		*c_i = a_i;
	} else {
		const double z2r   = 2.0 * a_r;
		const double z2i   = 2.0 * a_i;
		const double zd    = mc_cosh(z2r) + mc_cos(z2i);
		const double z2rsh = mc_sinh(z2r);
		if (mc_isinf(z2rsh) && mc_isinf(zd)) {
			*c_r = z2rsh > 0.0 ? 1.0 : -1.0;
			*c_i = z2i   > 0.0 ? 0.0 : -0.0;
		} else {
			*c_r = z2rsh / zd;
			*c_i = mc_sin(z2i) / zd;
		}
	}
#	endif
}

MC_TARGET_PROC void mc_ztanhl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
#	if MC_TARGET_EMBEDDED
	const long double z2r = 2.0L * a_r;
	const long double z2i = 2.0L * a_i;
	const long double zd  = mc_coshl(z2r) + mc_cosl(z2i);
	*c_r                  = mc_sinhl(z2r) / zd;
	*c_i                  = mc_sinl(z2i) / zd;
#	else
	if (mc_isinf(a_r)) {
		if (!mc_isfinite(a_i)) {
			*c_r = 1.0;
			*c_i = 0.0;
		} else {
			*c_r = 1.0;
			*c_i = mc_copysignl(0.0L, mc_sinl(2.0L * a_i));
		}
	} else if (mc_isnan(a_r) && a_i == 0.0L) {
		*c_r = a_r;
		*c_i = a_i;
	} else {
		const long double z2r   = 2.0L * a_r;
		const long double z2i   = 2.0L * a_i;
		const long double zd    = mc_coshl(z2r) + mc_cosl(z2i);
		const long double z2rsh = mc_sinhl(z2r);
		if (mc_isinf(z2rsh) && mc_isinf(zd)) {
			*c_r = z2rsh > 0.0L ? 1.0L : -1.0L;
			*c_i = z2i   > 0.0L ? 0.0L : -0.0L;
		} else {
			*c_r = z2rsh / zd;
			*c_i = mc_sinl(z2i) / zd;
		}
	}
#	endif
}

#endif /* !MC_ZTANH_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zatanh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_zadd.h>
#include <macadam/details/math/mc_zdiv.h>
#include <macadam/details/math/mc_zlog.h>
#include <macadam/details/math/mc_zmul.h>
#include <macadam/details/math/mc_zsub.h>

#ifndef MC_ZATANH_H
#define MC_ZATANH_H

#pragma mark - mc_zatanh -

MC_TARGET_PROC void mc_zatanhf(float * c_r, float * c_i
	, const float a_r, const float a_i
) {
	if (mc_isinf(a_i)) {
			*c_r = mc_copysignf(0.0f, a_r);
			*c_i = mc_copysignf(MCK_KF(MCK_PI_2), a_i);
	} else if (mc_isnan(a_i)) {
		if (mc_isinf(a_r) || a_r == 0.0f) {
			*c_r = mc_copysignf(0.0f, a_r);
			*c_i = a_i;
		} else {
			*c_r = a_i;
			*c_i = a_i;
		}
	} else if (mc_isnan(a_r)) {
		*c_r = a_r;
		*c_i = a_r;
	} else if (mc_isinf(a_r)) {
		*c_r = mc_copysignf(0.0f, a_r);
		*c_i = mc_copysignf(MCK_KF(MCK_PI_2), a_i);
	} else if (mc_fabsf(a_r) == 1.0f && a_i == 0.0f) {
		*c_r = mc_copysignf(MCK_INF, a_r);
		*c_i = mc_copysignf(0.0f, a_i);
	} else {
		float z1p_r, z1p_i;
		float z1m_r, z1m_i;

		mc_zaddf(&z1p_r, &z1p_i, 1.0f, 0.0f, a_r, a_i);
		mc_zsubf(&z1m_r, &z1m_i, 1.0f, 0.0f, a_r, a_i);

		mc_zdivf(c_r, c_i, z1p_r, z1p_i, z1m_r, z1m_i);
		mc_zlogf(c_r, c_i, *c_r, *c_i);
		mc_zmulf(c_r, c_i, *c_r, *c_i, 0.5f, 0.0f);

		*c_r = mc_copysignf(*c_r, a_r);
		*c_i = mc_copysignf(*c_i, a_i);
	}
}

MC_TARGET_PROC void mc_zatanh(double * c_r, double * c_i
	, const double a_r, const double a_i
) {
	if (mc_isinf(a_i)) {
			*c_r = mc_copysign(0.0, a_r);
			*c_i = mc_copysign(MCK_K(MCK_PI_2), a_i);
	} else if (mc_isnan(a_i)) {
		if (mc_isinf(a_r) || a_r == 0.0) {
			*c_r = mc_copysign(0.0, a_r);
			*c_i = a_i;
		} else {
			*c_r = a_i;
			*c_i = a_i;
		}
	} else if (mc_isnan(a_r)) {
		*c_r = a_r;
		*c_i = a_r;
	} else if (mc_isinf(a_r)) {
		*c_r = mc_copysign(0.0, a_r);
		*c_i = mc_copysign(MCK_K(MCK_PI_2), a_i);
	} else if (mc_fabs(a_r) == 1.0 && a_i == 0.0) {
		*c_r = mc_copysign(MCK_INF, a_r);
		*c_i = mc_copysign(0.0, a_i);
	} else {
		double z1p_r, z1p_i;
		double z1m_r, z1m_i;

		mc_zadd(&z1p_r, &z1p_i, 1.0, 0.0, a_r, a_i);
		mc_zsub(&z1m_r, &z1m_i, 1.0, 0.0, a_r, a_i);

		mc_zdiv(c_r, c_i, z1p_r, z1p_i, z1m_r, z1m_i);
		mc_zlog(c_r, c_i, *c_r, *c_i);
		mc_zmul(c_r, c_i, *c_r, *c_i, 0.5, 0.0);

		*c_r = mc_copysign(*c_r, a_r);
		*c_i = mc_copysign(*c_i, a_i);
	}
}

MC_TARGET_PROC void mc_zatanhl(long double * c_r, long double * c_i
	, const long double a_r, const long double a_i
) {
	if (mc_isinf(a_i)) {
			*c_r = mc_copysignl(0.0L, a_r);
			*c_i = mc_copysignl(MCK_KL(MCK_PI_2), a_i);
	} else if (mc_isnan(a_i)) {
		if (mc_isinf(a_r) || a_r == 0.0L) {
			*c_r = mc_copysignl(0.0L, a_r);
			*c_i = a_i;
		} else {
			*c_r = a_i;
			*c_i = a_i;
		}
	} else if (mc_isnan(a_r)) {
		*c_r = a_r;
		*c_i = a_r;
	} else if (mc_isinf(a_r)) {
		*c_r = mc_copysignl(0.0L, a_r);
		*c_i = mc_copysignl(MCK_KL(MCK_PI_2), a_i);
	} else if (mc_fabsl(a_r) == 1.0L && a_i == 0.0L) {
		*c_r = mc_copysignl(MCK_INF, a_r);
		*c_i = mc_copysignl(0.0L, a_i);
	} else {
		long double z1p_r, z1p_i;
		long double z1m_r, z1m_i;

		mc_zaddl(&z1p_r, &z1p_i, 1.0L, 0.0L, a_r, a_i);
		mc_zsubl(&z1m_r, &z1m_i, 1.0L, 0.0L, a_r, a_i);

		mc_zdivl(c_r, c_i, z1p_r, z1p_i, z1m_r, z1m_i);
		mc_zlogl(c_r, c_i, *c_r, *c_i);
		mc_zmull(c_r, c_i, *c_r, *c_i, 0.5L, 0.0L);

		*c_r = mc_copysignl(*c_r, a_r);
		*c_i = mc_copysignl(*c_i, a_i);
	}
}

#endif /* !MC_ZATANH_H */

/* EOF */
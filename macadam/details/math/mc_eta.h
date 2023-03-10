//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_eta.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_powm1.h>
#include <macadam/details/math/mc_riemann_zeta.h>

#ifndef MC_ETA_H
#define MC_ETA_H

#pragma mark - mc_eta -

MC_TARGET_FUNC float mc_etaf(const float x)
{
	if (mc_isnan(x)) {
		return MCK_NAN;
	}
	if (x == 1.0f) {
		return MCK_KF(MCK_LOGE2);
	}
	const float y = 1.0f - x;
	if (mc_fabsf(y) <= 2.000000000000000E-06f) {
		return MCK_KF(MCK_LOGE2) - (3.268629627000000E-02f * y + 1.598689037424310E-01f) * y;
	}
	const float zeta  = mc_riemann_zetaf(x);
	if (zeta == 0.0f || zeta == 1.0f) {
		return zeta;
	}
	if (mc_isnan(zeta)) {
		return MCK_NAN;
	}
	if (mc_isinf(zeta)) {
		return MCK_INFP;
	}
	const float powm1 = mc_powm1f(2.0f, 1.0f - x);
	if (mc_isnan(powm1)) {
		return MCK_NAN;
	}
	if (mc_isinf(powm1)) {
		return MCK_INFP;
	}
	return -powm1 * zeta;
}

MC_TARGET_FUNC double mc_eta(const double x)
{
	if (mc_isnan(x)) {
		return MCK_NAN;
	}
	if (x == 1.0) {
		return MCK_K(MCK_LOGE2);
	}
	const double y = 1.0 - x;
	if (mc_fabs(y) <= 2.000000000000000E-06) {
		return MCK_K(MCK_LOGE2) - (3.268629627000000E-02 * y + 1.598689037424310E-01) * y;
	}
	const double zeta  = mc_riemann_zeta(x);
	if (zeta == 0.0 || zeta == 1.0) {
		return zeta;
	}
	if (mc_isnan(zeta)) {
		return MCK_NAN;
	}
	if (mc_isinf(zeta)) {
		return MCK_INFP;
	}
	const double powm1 = mc_powm1(2.0, 1.0 - x);
	if (mc_isnan(powm1)) {
		return MCK_NAN;
	}
	if (mc_isinf(powm1)) {
		return MCK_INFP;
	}
	return -powm1 * zeta;
}

MC_TARGET_FUNC long double mc_etal(const long double x)
{
	if (mc_isnan(x)) {
		return MCK_NAN;
	}
	if (x == 1.0L) {
		return MCK_KL(MCK_LOGE2);
	}
	const long double y = 1.0L - x;
	if (mc_fabsl(y) <= 2.000000000000000E-06L) {
		return MCK_KL(MCK_LOGE2) - (3.268629627000000E-02L * y + 1.598689037424310E-01L) * y;
	}
	const long double zeta  = mc_riemann_zetal(x);
	if (zeta == 0.0L || zeta == 1.0L) {
		return zeta;
	}
	if (mc_isnan(zeta)) {
		return MCK_NAN;
	}
	if (mc_isinf(zeta)) {
		return MCK_INFP;
	}
	const long double powm1 = mc_powm1l(2.0L, 1.0L - x);
	if (mc_isnan(powm1)) {
		return MCK_NAN;
	}
	if (mc_isinf(powm1)) {
		return MCK_INFP;
	}
	return -powm1 * zeta;
}

#endif /* !MC_ETA_H */

/* EOF */
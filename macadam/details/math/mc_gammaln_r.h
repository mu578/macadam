//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_gammaln_r.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_gammaln.h>

#ifndef MC_GAMMALN_R_H
#define MC_GAMMALN_R_H

#pragma mark - mc_gammaln -

MC_TARGET_PROC float mc_gammalnf_r(const float x, int * psigngam)
{
//!# Computes log(|gamma(x)|).
	if (mc_isnan(x)) {
		return x;
	}
	if (mc_isinf(x)) {
		return x;
	}
	if (x <= 0.0f && mc_fisintf(x)) {
		return MCK_INFP;
	}
	if (x == 0.0f) {
		return MCK_INFP;
	}
	if (x == 1.0f || x == 2.0f) {
		return 0.0f;
	}
	const float y = mc_fabsf(x);
	if (y <= FLT_MIN) {
		//!# @todo gamma sign.
		return -mc_logf(y);
	}
	if (y > (MCLIMITS_MAXF / mc_logf(MCLIMITS_MAXF))) {
		return MCK_INFP;
	}
	return mc_gammalnf_approx0(x, psigngam);
}

MC_TARGET_PROC double mc_gammaln_r(const double x, int * psigngam)
{
//!# Computes log(|gamma(x)|).
	if (mc_isnan(x)) {
		return x;
	}
	if (mc_isinf(x)) {
		return x;
	}
	if (x <= 0.0 && mc_fisint(x)) {
		return MCK_INFP;
	}
	if (x == 0.0) {
		return MCK_INFP;
	}
	if (x == 1.0 || x == 2.0) {
		return 0.0;
	}
	const double y = mc_fabs(x);
	if (y <= DBL_MIN) {
		//!# @todo gamma sign.
		return -mc_log(y);
	}
	if (y > (MCLIMITS_MAX / mc_log(MCLIMITS_MAX))) {
		return MCK_INFP;
	}
	return mc_gammaln_approx0(x, psigngam);
}

MC_TARGET_PROC long double mc_gammalnl_r(const long double x, int * psigngam)
{
//!# Computes log(|gamma(x)|).
	if (mc_isnan(x)) {
		return x;
	}
	if (mc_isinf(x)) {
		return x;
	}
	if (x <= 0.0L && mc_fisintl(x)) {
		return MCK_INFP;
	}
	if (x == 0.0L) {
		return MCK_INFP;
	}
	if (x == 1.0L || x == 2.0L) {
		return 0.0L;
	}
	const long double y = mc_fabsl(x);
	if (y <= LDBL_MIN) {
		//!# @todo gamma sign.
		return -mc_logl(y);
	}
	if (y > (MCLIMITS_MAXL / mc_logl(MCLIMITS_MAXL))) {
		return MCK_INFP;
	}
	return mc_gammalnl_approx0(x, psigngam);
}

#endif /* !MC_GAMMALN_R_H */

/* EOF */
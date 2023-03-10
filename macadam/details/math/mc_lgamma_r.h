//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lgamma_r.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_gammaln.h>

#ifndef MC_LGAMMA_R_H
#define MC_LGAMMA_R_H

#pragma mark - mc_lgamma_r -

MC_TARGET_PROC float mc_lgammaf_r(const float x, int * psigngam)
{
//!# Computes log(|gamma(x)|).
#	if MC_TARGET_EMBEDDED || MC_TARGET_MSVC_CPP
	float g        = mc_gammalnf_r(x, psigngam);
	mc_gammasign_s = *psigngam;
	return g;
#	else
#	if MC_TARGET_CPP98
	float g        = ::lgammaf_r(x, psigngam);
	mc_gammasign_s = signgam;
	return g;
#	else
	float g        = lgammaf_r(x, psigngam);
	mc_gammasign_s = signgam;
	return g;
#	endif
#	endif
}

MC_TARGET_PROC double mc_lgamma_r(const double x, int * psigngam)
{
//!# Computes log(|gamma(x)|).
#	if MC_TARGET_EMBEDDED || MC_TARGET_MSVC_CPP
	double g       = mc_gammaln_r(x, psigngam);
	mc_gammasign_s = *psigngam;
	return g;
#	else
#	if MC_TARGET_CPP98
	double g       = ::lgamma_r(x, psigngam);
	mc_gammasign_s = signgam;
	return g;
#	else
	double g       = lgamma_r(x, psigngam);
	mc_gammasign_s = signgam;
	return g;
#	endif
#	endif
}

MC_TARGET_PROC long double mc_lgammal_r(const long double x, int * psigngam)
{
//!# Computes log(|gamma(x)|).
#	if MC_TARGET_EMBEDDED || MC_TARGET_MSVC_CPP
	long double g  = mc_gammalnl_r(x, psigngam);
	mc_gammasign_s = *psigngam;
	return g;
#	else
#	if MC_TARGET_CPP98
	long double g  = ::lgammal_r(x, psigngam);
	mc_gammasign_s = signgam;
	return g;
#	else
	long double g  = lgammal_r(x, psigngam);
	mc_gammasign_s = signgam;
	return g;
#	endif
#	endif
}

#endif /* !MC_LGAMMA_R_H */

/* EOF */
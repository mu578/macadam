//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lgamma.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_gammaln.h>

#ifndef MC_LGAMMA_H
#define MC_LGAMMA_H

#pragma mark - mc_lgamma -

MC_TARGET_FUNC float mc_lgammaf(const float x)
{
//!# Computes log(|gamma(x)|).
#	if MC_TARGET_EMBEDDED || MC_TARGET_MSVC_CPP
	return mc_gammalnf(x);
#	else
//!# IEEE Std 1003.1, 2004 lgamma implementation arbitrary
//!# trigger exceptions instead of mathematical rationality.
	if (x == 0.0f || (x < 0.0f && mc_fisintf(x))) {
		return MCK_INFP;
	}
#	if MC_TARGET_CPP98
	float g        = ::lgammaf(x);
	mc_gammasign_s = signgam;
	return g;
#	else
	float g        = lgammaf(x);
	mc_gammasign_s = signgam;
	return g;
#	endif
#	endif
}

MC_TARGET_FUNC double mc_lgamma(const double x)
{
//!# Computes log(|gamma(x)|).
#	if MC_TARGET_EMBEDDED || MC_TARGET_MSVC_CPP
	return mc_gammaln(x);
#	else
//!# IEEE Std 1003.1, 2004 lgamma implementation arbitrary
//!# trigger exceptions instead of mathematical rationality.
	if (x == 0.0 || (x < 0.0 && mc_fisint(x))) {
		return MCK_INFP;
	}
#	if MC_TARGET_CPP98
	double g       = ::lgamma(x);
	mc_gammasign_s = signgam;
	return g;
#	else
	double g       = lgamma(x);
	mc_gammasign_s = signgam;
	return g;
#	endif
#	endif
}

MC_TARGET_FUNC long double mc_lgammal(const long double x)
{
//!# Computes log(|gamma(x)|).
#	if MC_TARGET_EMBEDDED || MC_TARGET_MSVC_CPP
	return mc_gammalnl(x);
#	else
//!# IEEE Std 1003.1, 2004 lgamma implementation arbitrary
//!# trigger exceptions instead of mathematical rationality.
	if (x == 0.0L || (x < 0.0L && mc_fisintl(x))) {
		return MCK_INFP;
	}
#	if MC_TARGET_CPP98
	long double g  = ::lgammal(x);
	mc_gammasign_s = signgam;
	return g;
#	else
	long double g  = lgammal(x);
	mc_gammasign_s = signgam;
	return g;
#	endif
#	endif
}

#endif /* !MC_LGAMMA_H */

/* EOF */
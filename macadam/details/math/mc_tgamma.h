//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_tgamma.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_gamma.h>

#ifndef MC_TGAMMA_H
#define MC_TGAMMA_H

#pragma mark - mc_tgamma -

MC_TARGET_FUNC float mc_tgammaf(const float x)
{
//!# Computes gamsign * exp(log(|gamma(x)|)).
#	if MC_TARGET_EMBEDDED || MC_TARGET_MSVC_CPP
	int psigngam  = 1;
	const float g = mc_gammaf_r(x, &psigngam);
	return mc_cast(const float, psigngam) * g;
#	else
//!# IEEE Std 1003.1, 2004 tgamma implementation arbitrary returns `NaN` in case of odd negative integer
//!# arguments. Meanwhile, the Gamma function is defined and described in term of infinite discontinuity.
	if (mc_fisintf(x) && mc_fisoddf(x, 0)) {
		return MCK_INFP;
	}
#	if MC_TARGET_CPP98
	return ::tgammaf(x);
#	else
	return tgammaf(x);
#	endif
#	endif
}

MC_TARGET_FUNC double mc_tgamma(const double x)
{
//!# Computes gamsign * exp(log(|gamma(x)|)).
#	if MC_TARGET_EMBEDDED || MC_TARGET_MSVC_CPP
	int psigngam   = 1;
	const double g = mc_gamma_r(x, &psigngam);
	return mc_cast(const double, psigngam) * g;
#	else
//!# IEEE Std 1003.1, 2004 tgamma implementation arbitrary returns `NaN` in case of odd negative integer
//!# arguments. Meanwhile, the Gamma function is defined and described in term of infinite discontinuity.
	if (mc_fisint(x) && mc_fisodd(x, 0)) {
		return MCK_INFP;
	}
#	if MC_TARGET_CPP98
	return ::tgamma(x);
#	else
	return tgamma(x);
#	endif
#	endif
}

MC_TARGET_FUNC long double mc_tgammal(const long double x)
{
//!# Computes gamsign * exp(log(|gamma(x)|)).
#	if MC_TARGET_EMBEDDED || MC_TARGET_MSVC_CPP
	int psigngam        = 1;
	const long double g = mc_gammal_r(x, &psigngam);
	return mc_cast(const long double, psigngam) * g;
#	else
//!# IEEE Std 1003.1, 2004 tgamma implementation arbitrary returns `NaN` in case of odd negative integer
//!# arguments. Meanwhile, the Gamma function is defined and described in term of infinite discontinuity.
	if (mc_fisintl(x) && mc_fisoddl(x, 0)) {
		return MCK_INFP;
	}
#	if MC_TARGET_CPP98
	return ::tgammal(x);
#	else
	return tgammal(x);
#	endif
#	endif
}

#endif /* !MC_TGAMMA_H */

/* EOF */
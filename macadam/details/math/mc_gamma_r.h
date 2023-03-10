//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_gamma_r.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_gamma.h>

#ifndef MC_GAMMA_R_H
#define MC_GAMMA_R_H

#pragma mark - mc_gamma -

MC_TARGET_PROC float mc_gammaf_r(const float x, int * psigngam)
{
//!# Computes exp(log(|gamma(x)|)).
	return mc_gammaf_approx1(x, psigngam);
}

MC_TARGET_PROC double mc_gamma_r(const double x, int * psigngam)
{
//!# Computes exp(log(|gamma(x)|)).
	return mc_gamma_approx1(x, psigngam);
}

MC_TARGET_PROC long double mc_gammal_r(const long double x, int * psigngam)
{
//!# Computes exp(log(|gamma(x)|)).
	return mc_gammal_approx1(x, psigngam);
}

#endif /* !MC_GAMMA_R_H */

/* EOF */
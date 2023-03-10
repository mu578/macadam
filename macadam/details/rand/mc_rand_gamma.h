//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rand_gamma.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/rand/mc_rand_lgamma.h>

#ifndef MC_RAND_GAMMA_H
#define MC_RAND_GAMMA_H

#pragma mark - mc_rand_gamma -

MC_TARGET_FUNC float mc_rand_gammaf(const float a, const float l)
{
//!# Gamma(alpha,lambda) generator using Marsaglia and Tsang method.
//!# Gamma RNG a=alpha=shape, l=lambda=scale.
	return mc_expf(mc_rand_lgammaf(a, l));
}

MC_TARGET_FUNC double mc_rand_gammaff(const float a, const float l)
{
//!# Gamma(alpha,lambda) generator using Marsaglia and Tsang method.
//!# Gamma RNG a=alpha=shape, l=lambda=scale.
	return mc_exp(mc_rand_lgammaff(a, l));
}

MC_TARGET_FUNC double mc_rand_gamma(const double a, const double l)
{
//!# Gamma(alpha,lambda) generator using Marsaglia and Tsang method.
//!# Gamma RNG a=alpha=shape, l=lambda=scale.
	return mc_exp(mc_rand_lgamma(a, l));
}

MC_TARGET_FUNC long double mc_rand_gammal(const long double a, const long double l)
{
//!# Gamma(alpha,lambda) generator using Marsaglia and Tsang method.
//!# Gamma RNG a=alpha=shape, l=lambda=scale.
	return mc_expl(mc_rand_lgammal(a, l));
}

#endif /* !MC_RAND_GAMMA_H */

/* EOF */
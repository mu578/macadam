//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rand_chisquare.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/rand/mc_rand_gamma.h>

#ifndef MC_RAND_CHISQUARE_H
#define MC_RAND_CHISQUARE_H

#pragma mark - mc_rand_chisquare -

MC_TARGET_FUNC float mc_rand_chisquaref(const float df)
{
//!# Chisquare distribution generator. df=degree of freedom.
	return 2.0f * mc_rand_gammaf(df * 0.5f, 1.0f);
}

MC_TARGET_FUNC double mc_rand_chisquareff(const float df)
{
//!# Chisquare distribution generator. df=degree of freedom.
	return 2.0 * mc_rand_gammaff(df * 0.5f, 1.0f);
}

MC_TARGET_FUNC double mc_rand_chisquare(const double df)
{
//!# Chisquare distribution generator. df=degree of freedom.
	return 2.0 * mc_rand_gamma(df * 0.5, 1.0);
}

MC_TARGET_FUNC long double mc_rand_chisquarel(const long double df)
{
//!# Chisquare distribution generator. df=degree of freedom.
	return 2.0L * mc_rand_gammal(df * 0.5L, 1.0L);
}

#endif /* !MC_RAND_CHISQUARE_H */

/* EOF */
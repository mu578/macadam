//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_gamma_p.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_igamma.h>

#ifndef MC_GAMMAP_H
#define MC_GAMMAP_H

#pragma mark - mc_gamma_p -

MC_TARGET_FUNC float mc_gamma_pf(const float a, const float z)
{
	return mc_igamma_pf_approx2(a, z);
}

MC_TARGET_FUNC double mc_gamma_p(const double a, const double z)
{
	return mc_igamma_p_approx2(a, z);
}

MC_TARGET_FUNC long double mc_gamma_pl(long double a, long double z)
{
	return mc_igamma_pl_approx2(a, z);
}

#endif /* !MC_GAMMAP_H */

/* EOF */
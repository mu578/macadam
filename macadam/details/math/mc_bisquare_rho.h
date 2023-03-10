//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_bisquare_rho.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_raise3.h>

#ifndef MC_BISQUARE_RHO_H
#define MC_BISQUARE_RHO_H

#pragma mark - mc_bisquare_rhoc -

MC_TARGET_FUNC float mc_bisquare_rhocf(const float r, const float c)
{
	if (mc_fabsf(r) <= c) {
		return (mc_raise2f(c) / 6.0f) * mc_raise3f(1.0f - (1.0f - mc_raise2f(r / c)));
	}
	return mc_raise2f(c) / 6.0f;
}

MC_TARGET_FUNC double mc_bisquare_rhoc(const double r, const double c)
{
	if (mc_fabs(r) <= c) {
		return (mc_raise2(c) / 6.0) * mc_raise3(1.0 - (1.0 - mc_raise2(r / c)));
	}
	return mc_raise2(c) / 6.0;
}

MC_TARGET_FUNC long double mc_bisquare_rhocl(const long double r, const long double c)
{
	if (mc_fabsl(r) <= c) {
		return (mc_raise2l(c) / 6.0L) * mc_raise3l(1.0L - (1.0L - mc_raise2l(r / c)));
	}
	return mc_raise2l(c) / 6.0L;
}

#pragma mark - mc_bisquare_rho -

MC_TARGET_FUNC float mc_bisquare_rhof(const float r)
{
	const float c = 4.685f;
	return mc_bisquare_rhocf(r, c);
}

MC_TARGET_FUNC double mc_bisquare_rho(const double r)
{
	const double c = 4.685;
	return mc_bisquare_rhoc(r, c);
}

MC_TARGET_FUNC long double mc_bisquare_rhol(const long double r)
{
	const long double c = 4.685L;
	return mc_bisquare_rhocl(r, c);
}

#endif /* !MC_BISQUARE_RHO_H */

/* EOF */
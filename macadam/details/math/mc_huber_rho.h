//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_huber_rho.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_HUBER_RHO_H
#define MC_HUBER_RHO_H

#pragma mark - mc_huber_rhoc -

MC_TARGET_FUNC float mc_huber_rhocf(const float r, const float c)
{
	if (mc_fabsf(r) <= c) {
		return 0.5f * mc_raise2f(r);
	}
	return c * mc_fabsf(r) - 0.5f * mc_raise2f(c);
}

MC_TARGET_FUNC double mc_huber_rhoc(const double r, const double c)
{
	if (mc_fabs(r) <= c) {
		return 0.5 * mc_raise2(r);
	}
	return c * mc_fabs(r) - 0.5 * mc_raise2(c);
}

MC_TARGET_FUNC long double mc_huber_rhocl(const long double r, const long double c)
{
	if (mc_fabsl(r) <= c) {
		return 0.5L * mc_raise2l(r);
	}
	return c * mc_fabsl(r) - 0.5L * mc_raise2l(c);
}

#pragma mark - mc_huber_rho -

MC_TARGET_FUNC float mc_huber_rhof(const float r)
{
	const float c = 1.345f;
	return mc_huber_rhocf(r, c);
}

MC_TARGET_FUNC double mc_huber_rho(const double r)
{
	const double c = 1.345;
	return mc_huber_rhoc(r, c);
}

MC_TARGET_FUNC long double mc_huber_rhol(const long double r)
{
	const long double c = 1.345L;
	return mc_huber_rhocl(r, c);
}

#endif /* !MC_HUBER_RHO_H */

/* EOF */
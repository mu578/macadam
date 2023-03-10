//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_bisquare_psi.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_BISQUARE_PSI_H
#define MC_BISQUARE_PSI_H

#pragma mark - mc_bisquare_psic -

MC_TARGET_FUNC float mc_bisquare_psicf(const float r, const float c)
{
	if (mc_fabsf(r) <= c) {
		return r * mc_raise2f(mc_raise2f(1.0f - (r / c)));
	}
	return 0.0f;
}

MC_TARGET_FUNC double mc_bisquare_psic(const double r, const double c)
{
	if (mc_fabs(r) <= c) {
		return r * mc_raise2(mc_raise2(1.0 - (r / c)));
	}
	return 0.0;
}

MC_TARGET_FUNC long double mc_bisquare_psicl(const long double r, const long double c)
{
	if (mc_fabsl(r) <= c) {
		return r * mc_raise2l(mc_raise2l(1.0L - (r / c)));
	}
	return 0.0L;
}

#pragma mark - mc_bisquare_psi -

MC_TARGET_FUNC float mc_bisquare_psif(const float r)
{
	const float c = 4.685f;
	return mc_bisquare_psicf(r, c);
}

MC_TARGET_FUNC double mc_bisquare_psi(const double r)
{
	const double c = 4.685;
	return mc_bisquare_psic(r, c);
}

MC_TARGET_FUNC long double mc_bisquare_psil(const long double r)
{
	const long double c = 4.685L;
	return mc_bisquare_psicl(r, c);
}

#endif /* !MC_BISQUARE_PSI_H */

/* EOF */
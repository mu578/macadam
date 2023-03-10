//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fmodpi.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fisnear.h>
#include <macadam/details/math/mc_fmod.h>

#ifndef MC_FMODPI_H
#define MC_FMODPI_H

#pragma mark - mc_fmodpi -

MC_TARGET_FUNC float mc_fmodpif(const float x)
{
//!# [-pi, +pi]. x mod pi.
	const float p1 = MCK_KF(MCK_PI);
	const float p2 = MCK_KF(MCK_2PI);
	const float y  = (x + p1) / p2;
	return (x - p2 * mc_floorf(y)) + 0.0f;
}

MC_TARGET_FUNC double mc_fmodpi(const double x)
{
//!# [-pi, +pi]. x mod pi.
	const double p1 = MCK_K(MCK_PI);
	const double p2 = MCK_K(MCK_2PI);
	const double y  = (x + p1) / p2;
	return (x - p2 * mc_floor(y)) + 0.0;
}

MC_TARGET_FUNC long double mc_fmodpil(const long double x)
{
//!# [-pi, +pi]. x mod pi.
	const long double p1 = MCK_KL(MCK_PI);
	const long double p2 = MCK_KL(MCK_2PI);
	const long double y  = (x + p1) / p2;
	return (x - p2 * mc_floorl(y)) + 0.0L;
}

#endif /* !MC_FMODPI_H */

/* EOF */
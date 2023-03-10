//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fmod2pi.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fisnear.h>
#include <macadam/details/math/mc_fmod.h>

#ifndef MC_FMOD2PI_H
#define MC_FMOD2PI_H

#pragma mark - mc_fmod2pi -

MC_TARGET_FUNC float mc_fmod2pif(const float x)
{
//!# [+0, +2pi].
	const float p = MCK_KF(MCK_2PI);
	if (x == 0.0f || mc_fisnearf(x, p, 2)) {
		return 0.0f;
	}
	// err' ~= +(30.0f * MCLIMITS_EPSILONF);
	const float m = mc_fmodf(x, p);
	return m < 0.0f ? m + p : m + 0.0f;
}

MC_TARGET_FUNC double mc_fmod2pi(const double x)
{
//!# [+0, +2pi].
	const double p = MCK_K(MCK_2PI);
	if (x == 0.0 || mc_fisnear(x, p, 2)) {
		return 0.0;
	}
	// err' ~= -(22.0 * MCLIMITS_EPSILON);
	const double m = mc_fmod(x, p);
	return m < 0.0 ? m + p : m + 0.0;
}

MC_TARGET_FUNC long double mc_fmod2pil(const long double x)
{
//!# [+0, +2pi].
	const long double p = MCK_KL(MCK_2PI);
	if (x == 0.0L || mc_fisnearl(x, p, 2)) {
		return 0.0L;
	}
	const long double m = mc_fmodl(x, p);
	return m < 0.0L ? m + p : m + 0.0L;
}

#endif /* !MC_FMOD2PI_H */

/* EOF */
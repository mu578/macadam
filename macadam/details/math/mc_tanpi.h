//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_tanpi.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_fmod.h>
#include <macadam/details/math/mc_tan.h>

#ifndef MC_TANPI_H
#define MC_TANPI_H

#pragma mark - mc_tanpi -

MC_TARGET_FUNC float mc_tanpif(const float x)
{
	float y = x;
	y       = mc_fmodf(y, 1.0f);
	y       = y <= -0.5f ? y + 1.0f : y > 0.5f ? y - 1.0f : y;
	if (y == 1.0f) {
		return mc_copysignf(0.0f, y);
	}
	if (y == 0.5f) {
		return mc_copysignf(MCK_INF, y);
	}
	const float piy = MCK_KF(MCK_PI) * y;
	return mc_tanf(piy);
}

MC_TARGET_FUNC double mc_tanpi(const double x)
{
	double y = x;
	y        = mc_fmod(y, 1.0);
	y        = y <= -0.5 ? y + 1.0 : y > 0.5 ? y - 1.0 : y;
	if (y == 1.0) {
		return mc_copysign(0.0, y);
	}
	if (y == 0.5) {
		return mc_copysign(MCK_INF, y);
	}
	const double piy = MCK_K(MCK_PI) * y;
	return mc_tan(piy);
}

MC_TARGET_FUNC long double mc_tanpil(const long double x)
{
	long double y = x;
	y             = mc_fmodl(y, 1.0L);
	y             = y <= -0.5L ? y + 1.0L : y > 0.5L ? y - 1.0L : y;
	if (y == 1.0L) {
		return mc_copysignl(0.0L, y);
	}
	if (y == 0.5L) {
		return mc_copysignl(MCK_INF, y);
	}
	const long double piy = MCK_KL(MCK_PI) * y;
	return mc_tanl(piy);
}

#endif /* !MC_TANPI_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lerp.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_isnan.h>

#ifndef MC_LERP_H
#define MC_LERP_H

#pragma mark - mc_lerp -

MC_TARGET_FUNC float mc_lerpf(const float x, const float y, const float z)
{
	if (mc_isnan(x) || mc_isnan(y) || mc_isnan(z)) {
		return MCK_NAN;
	}
	if ((x <= 0.0f && y >= 0.0f) || (x >= 0.0f && y <= 0.0f)) {
		return (x * (1.0f - z)) + (y * z);
	} else if (z == 1.0f) {
		return y;
	}
	const float w = x + z * (y - x);
	if (z > 1.0f && y > x) {
		return y < w ? w : y;
	}
	return w < y ? w : y;
}

MC_TARGET_FUNC double mc_lerp(const double x, const double y, const double z)
{
	if (mc_isnan(x) || mc_isnan(y) || mc_isnan(z)) {
		return MCK_NAN;
	}
	if ((x <= 0.0 && y >= 0.0) || (x >= 0.0 && y <= 0.0)) {
		return (x * (1.0 - z)) + (y * z);
	} else if (z == 1.0) {
		return y;
	}
	const double w = x + z * (y - x);
	if (z > 1.0 && y > x) {
		return y < w ? w : y;
	}
	return w < y ? w : y;
}

MC_TARGET_FUNC long double mc_lerpl(const long double x, const long double y, const long double z)
{
	if (mc_isnan(x) || mc_isnan(y) || mc_isnan(z)) {
		return MCK_NAN;
	}
	if ((x <= 0.0L && y >= 0.0L) || (x >= 0.0L && y <= 0.0L)) {
		return (x * (1.0L - z)) + (y * z);
	} else if (z == 1.0L) {
		return y;
	}
	const long double w = x + z * (y - x);
	if (z > 1.0L && y > x) {
		return y < w ? w : y;
	}
	return w < y ? w : y;
}

#endif /* !MC_LERP_H */

/* EOF */
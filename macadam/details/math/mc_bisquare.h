//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_bisquare.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_BISQUARE_H
#define MC_BISQUARE_H

#pragma mark - mc_bisquare -

MC_TARGET_FUNC float mc_bisquaref(const float r, const float c, const float s)
{
//!# Tukey's bisquare function.
//!# Default settings c=4.685 and s=1 (scale)
	const float w = s != 0.0f ? 1.0f / s : 1.0f;
	const float h = mc_fabsf(r * w);
	if (h > c) {
		return 0.0f;
	}
	return mc_raise2f(1.0f - mc_raise2f(h / c));
}

MC_TARGET_FUNC double mc_bisquare(const double r, const double c, const double s)
{
//!# Tukey's bisquare function.
//!# Default settings c=4.685 and s=1 (scale)
	const double w = s != 0.0 ? 1.0 / s : 1.0;
	const double h = mc_fabs(r * w);
	if (h > c) {
		return 0.0;
	}
	return mc_raise2(1.0 - mc_raise2(h / c));
}

MC_TARGET_FUNC long double mc_bisquarel(const long double r, const long double c, const long double s)
{
//!# Tukey's bisquare function.
//!# Default settings c=4.685 and s=1 (scale)
	const long double w = s != 0.0L ? 1.0L / s : 1.0L;
	const long double h = mc_fabsl(r * w);
	if (h > c) {
		return 0.0L;
	}
	return mc_raise2l(1.0L - mc_raise2l(h / c));
}

#endif /* !MC_BISQUARE_H */

/* EOF */
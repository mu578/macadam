//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sqrt1pm1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_SQRT1PM1_H
#define MC_SQRT1PM1_H

#pragma mark - mc_sqrt1pm1 -

MC_TARGET_FUNC float mc_sqrt1pm1f(const float x)
{
	if (mc_fabsf(x) > 0.75f) {
		return mc_sqrtf(1.0f + x) - 1.0f;
	}
	return x / (1.0f + mc_sqrtf(1.0f + x));
}

MC_TARGET_FUNC double mc_sqrt1pm1(const double x)
{
	if (mc_fabs(x) > 0.75) {
		return mc_sqrt(1.0 + x) - 1.0;
	}
	return x / (1.0 + mc_sqrt(1.0 + x));
}

MC_TARGET_FUNC long double mc_sqrt1pm1l(const long double x)
{
	if (mc_fabsl(x) > 0.75L) {
		return mc_sqrtl(1.0L + x) - 1.0L;
	}
	return x / (1.0L + mc_sqrtl(1.0L + x));
}

#endif /* !MC_SQRT1PM1_H */

/* EOF */
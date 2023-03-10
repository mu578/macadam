//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_hypot2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_HYPOT2_H
#define MC_HYPOT2_H

#pragma mark - mc_hypot2 -

MC_TARGET_FUNC float mc_hypot2f(const float x, const float y)
{
	float a = mc_fabsf(x);
	float b = mc_fabsf(y);
	float c = 0.0f;
	if (a > b) {
		c = b / a;
		c = a * mc_sqrtf(1.0f + c * c);
	} else if (b > 0.0f) {
		c = a / b;
		c = b * mc_sqrtf(1.0f + c * c);
	}
	return c;
}

MC_TARGET_FUNC double mc_hypot2(const double x, const double y)
{
	double a = mc_fabs(x);
	double b = mc_fabs(y);
	double c = 0.0;
	if (a > b) {
		c = b / a;
		c = a * mc_sqrt(1.0 + c * c);
	} else if (b > 0.0) {
		c = a / b;
		c = b * mc_sqrt(1.0 + c * c);
	}
	return c;
}

MC_TARGET_FUNC long double mc_hypot2l(const long double x, const long double y)
{
	long double a = mc_fabsl(x);
	long double b = mc_fabsl(y);
	long double c = 0.0L;
	if (a > b) {
		c = b / a;
		c = a * mc_sqrtl(1.0L + c * c);
	} else if (b > 0.0L) {
		c = a / b;
		c = b * mc_sqrtl(1.0L + c * c);
	}
	return c;
}

#endif /* !MC_HYPOT2_H */

/* EOF */
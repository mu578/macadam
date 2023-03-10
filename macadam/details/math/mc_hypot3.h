//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_hypot3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_HYPOT3_H
#define MC_HYPOT3_H

#pragma mark - mc_hypot3 -

MC_TARGET_FUNC float mc_hypot3f(const float x, const float y, const float z)
{
	float a = mc_fabsf(x);
	float b = mc_fabsf(y);
	float c = mc_fabsf(z);
	float d = a < b ? b < c ? c : b : a < c ? c : a;
	return (d != 0.0f
		? (d * mc_sqrtf((a / d) * (a / d) + (b / d) * (b / d) + (c / d) * (c / d)))
		: 0.0f
	);
}

MC_TARGET_FUNC double mc_hypot3(const double x, const double y, const double z)
{
	double a = mc_fabs(x);
	double b = mc_fabs(y);
	double c = mc_fabs(z);
	double d = a < b ? b < c ? c : b : a < c ? c : a;
	return (d != 0.0
		? (d * mc_sqrt((a / d) * (a / d) + (b / d) * (b / d) + (c / d) * (c / d)))
		: 0.0
	);
}

MC_TARGET_FUNC long double mc_hypot3l(const long double x, const long double y, const long double z)
{
	long double a = mc_fabsl(x);
	long double b = mc_fabsl(y);
	long double c = mc_fabsl(z);
	long double d = a < b ? b < c ? c : b : a < c ? c : a;
	return (d != 0.0L
		? (d * mc_sqrtl((a / d) * (a / d) + (b / d) * (b / d) + (c / d) * (c / d)))
		: 0.0L
	);
}

#endif /* !MC_HYPOT3_H */

/* EOF */
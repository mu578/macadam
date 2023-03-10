//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_exp10m1m1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_exp2m1.h>
#include <macadam/details/math/mc_fabs.h>

#ifndef MC_EXP10M1_H
#define MC_EXP10M1_H

#pragma mark - mc_exp10m1 -

MC_TARGET_FUNC float mc_exp10m1f(const float x)
{
	if (x == 0.0f) {
		return x;
	}
	if (mc_fabsf(x) > 0.3f) {
		return mc_exp10f(x) - 1.0f;
	}
	return mc_exp2m1f(MCK_KF(MCK_LOG210Q) * x + MCK_KF(MCK_LOG210R) * x);
}

MC_TARGET_FUNC double mc_exp10m1(const double x)
{
	if (x == 0.0) {
		return x;
	}
	if (mc_fabs(x) > 0.3) {
		return mc_exp10(x) - 1.0;
	}
	return mc_exp2m1(MCK_K(MCK_LOG210Q) * x + MCK_K(MCK_LOG210R) * x);
}

MC_TARGET_FUNC long double mc_exp10m1l(const long double x)
{
	if (x == 0.0L) {
		return x;
	}
	if (mc_fabsl(x) > 0.3L) {
		return mc_exp10l(x) - 1.0L;
	}
	return mc_exp2m1l(MCK_KL(MCK_LOG210Q) * x + MCK_KL(MCK_LOG210R) * x);
}

#endif /* !MC_EXP10M1_H */

/* EOF */
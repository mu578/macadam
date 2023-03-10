//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sigmoid.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_SIGMOID_H
#define MC_SIGMOID_H

#pragma mark - mc_sigmoid -

MC_TARGET_FUNC float mc_sigmoidf(const float x)
{
	const float x0 = +46.0f;
	const float x1 = -23.0f;
	if (mc_fabsf(x) > x0) {
		if (x > 0.0f) {
			return 1.0f;
		}
		return mc_expf(x);
	}
	if (x < x1) {
		const float y = mc_expf(x);
		return y - mc_raise2f(y);
	}
	return 1.0f / (mc_expf(-x) + 1.0f);
}

MC_TARGET_FUNC double mc_sigmoid(const double x)
{
	const double x0 = +46.0;
	const double x1 = -23.0;
	if (mc_fabs(x) > x0) {
		if (x > 0.0) {
			return 1.0;
		}
		return mc_exp(x);
	}
	if (x < x1) {
		const double y = mc_exp(x);
		return y - mc_raise2(y);
	}
	return 1.0 / (mc_exp(-x) + 1.0);
}

MC_TARGET_FUNC long double mc_sigmoidl(const long double x)
{
	const long double x0 = +46.0;
	const long double x1 = -23.0;
	if (mc_fabsl(x) > x0) {
		if (x > 0.0) {
			return 1.0;
		}
		return mc_expl(x);
	}
	if (x < x1) {
		const long double y = mc_expl(x);
		return y - mc_raise2l(y);
	}
	return 1.0L / (mc_expl(-x) + 1.0L);
}

#endif /* !MC_SIGMOID_H */

/* EOF */
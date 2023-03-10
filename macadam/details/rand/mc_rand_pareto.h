//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rand_pareto.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/rand/mc_rand_exponential.h>

#ifndef MC_RAND_PARETO_H
#define MC_RAND_PARETO_H

#pragma mark - mc_rand_pareto -

MC_TARGET_FUNC float mc_rand_paretof(const float a)
{
//!# Pareto distribution generator.
	float r = 0.0f;
	if (a != r) {
		const float x = mc_rand_exponentialf(1.0f) / a;
		r             = mc_expf(x) - 1.0f;
	}
	return r;
}

MC_TARGET_FUNC double mc_rand_paretoff(const float a)
{
//!# Pareto distribution generator.
	double r = 0.0;
	if (mc_cast(double, a) != r) {
		const double x = mc_rand_exponential(1.0) / mc_cast(const double, a);
		r              = mc_exp(x) - 1.0;
	}
	return r;
}

MC_TARGET_FUNC double mc_rand_pareto(const double a)
{
//!# Pareto distribution generator.
	double r = 0.0;
	if (a != r) {
		const double x = mc_rand_exponential(1.0) / a;
		r              = mc_exp(x) - 1.0;
	}
	return r;
}

MC_TARGET_FUNC long double mc_rand_paretol(const long double a)
{
//!# Pareto distribution generator.
	long double r = 0.0L;
	if (a != r) {
		const long double x = mc_rand_exponentiall(1.0L) / a;
		r                   = mc_expl(x) - 1.0L;
	}
	return r;
}

#endif /* !MC_RAND_PARETO_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_spvar1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/numa/mc_mssqr1xn.h>

#ifndef MC_SPVAR1XN_H
#define MC_SPVAR1XN_H

#pragma mark - mc_spvar1xn -

MC_TARGET_FUNC float mc_spvar1xnf(const int n, const float * x)
{
//!# Returns the sample variance of vector X.
	float mean, sumsq, scale;
	if (n < 2) {
		return 0.0f;
	}
	mc_mssqr1xnf(n, x, 0, &mean, &sumsq, &scale);
	return mc_raise2f(scale) * (sumsq / mc_cast_expr(const float, n - 1));
}

MC_TARGET_FUNC double mc_spvar1xnff(const int n, const float * x)
{
//!# Returns the sample variance of vector X.
	double mean, sumsq, scale;
	if (n < 2) {
		return 0.0;
	}
	mc_mssqr1xnff(n, x, 0, &mean, &sumsq, &scale);
	return mc_raise2(scale) * (sumsq / mc_cast_expr(const double, n - 1));
}

MC_TARGET_FUNC double mc_spvar1xn(const int n, const double * x)
{
//!# Returns the sample variance of vector X.
	double mean, sumsq, scale;
	if (n < 2) {
		return 0.0;
	}
	mc_mssqr1xn(n, x, 0, &mean, &sumsq, &scale);
	return mc_raise2(scale) * (sumsq / mc_cast_expr(const double, n - 1));
}

MC_TARGET_FUNC long double mc_spvar1xnl(const int n, const long double * x)
{
//!# Returns the sample variance of vector X.
	long double mean, sumsq, scale;
	if (n < 2) {
		return 0.0L;
	}
	mc_mssqr1xnl(n, x, 0, &mean, &sumsq, &scale);
	return mc_raise2l(scale) * (sumsq / mc_cast_expr(const long double, n - 1));
}

#endif /* !MC_SPVAR1XN_H */

/* EOF */
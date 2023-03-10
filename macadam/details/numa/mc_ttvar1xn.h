//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ttvar1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/numa/mc_mssqr1xn.h>

#ifndef MC_TTVAR1XN_H
#define MC_TTVAR1XN_H

#pragma mark - mc_ttvar1xn -

MC_TARGET_FUNC float mc_ttvar1xnf(const int n, const float * x)
{
//!# Returns the total variance of vector X.
	float mean, sumsq, scale;
	mc_mssqr1xnf(n, x, 0, &mean, &sumsq, &scale);
	return mc_raise2f(scale) * sumsq;
}

MC_TARGET_FUNC double mc_ttvar1xnff(const int n, const float * x)
{
//!# Returns the total variance of vector X.
	double mean, sumsq, scale;
	mc_mssqr1xnff(n, x, 0, &mean, &sumsq, &scale);
	return mc_raise2(scale) * sumsq;
}

MC_TARGET_FUNC double mc_ttvar1xn(const int n, const double * x)
{
//!# Returns the total variance of vector X.
	double mean, sumsq, scale;
	mc_mssqr1xn(n, x, 0, &mean, &sumsq, &scale);
	return mc_raise2(scale) * sumsq;
}

MC_TARGET_FUNC long double mc_ttvar1xnl(const int n, const long double * x)
{
//!# Returns the total variance of vector X.
	long double mean, sumsq, scale;
	mc_mssqr1xnl(n, x, 0, &mean, &sumsq, &scale);
	return mc_raise2l(scale) * sumsq;
}

#endif /* !MC_TTVAR1XN_H */

/* EOF */
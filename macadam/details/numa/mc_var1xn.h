//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_var1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/numa/mc_mssqr1xn.h>

#ifndef MC_VAR1XN_H
#define MC_VAR1XN_H

#pragma mark - mc_var1xn -

MC_TARGET_FUNC void mc_var1xnf(const int n, const float * x, const int b, float * var)
{
	float mean, sumsq, scale;

	*var = 0.0f;
	if (n > 1) {
		mc_mssqr1xnf(n, x, 0, &mean, &sumsq, &scale);
		*var = mc_raise2f(scale) * (sumsq / mc_cast(const float, (b ? (n - 1) : n)));
	}
}

MC_TARGET_FUNC void mc_var1xnff(const int n, const float * x, const int b, double * var)
{
	double mean, sumsq, scale;

	*var = 0.0;
	if (n > 1) {
		mc_mssqr1xnff(n, x, 0, &mean, &sumsq, &scale);
		*var = mc_raise2(scale) * (sumsq / mc_cast(const double, (b ? (n - 1) : n)));
	}
}

MC_TARGET_FUNC void mc_var1xn(const int n, const double * x, const int b, double * var)
{
	double mean, sumsq, scale;

	*var = 0.0;
	if (n > 1) {
		mc_mssqr1xn(n, x, 0, &mean, &sumsq, &scale);
		*var = mc_raise2(scale) * (sumsq / mc_cast(const double, (b ? (n - 1) : n)));
	}
}

MC_TARGET_FUNC void mc_var1xnl(const int n, const long double * x, const int b, long double * var)
{
	long double mean, sumsq, scale;

	*var = 0.0L;
	if (n > 1) {
		mc_mssqr1xnl(n, x, 0, &mean, &sumsq, &scale);
		*var = mc_raise2l(scale) * (sumsq / mc_cast(const long double, (b ? (n - 1) : n)));
	}
}

#endif /* !MC_VAR1XN_H */

/* EOF */
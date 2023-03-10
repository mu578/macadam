//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mstdd1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/numa/mc_mssqr1xn.h>

#ifndef MC_MSTDD1XN_H
#define MC_MSTDD1XN_H

#pragma mark - mc_mstdd1xn -

MC_TARGET_FUNC void mc_mstdd1xnf(const int n, const float * x, const int b, float * mean, float * stddev)
{
	float sumsq, scale;

	*mean   = 0.0f;
	*stddev = 0.0f;
	if (n > 1) {
		mc_mssqr1xnf(n, x, (b > 1 ? 1 : 0), mean, &sumsq, &scale);
		*stddev = scale * mc_sqrtf(sumsq / mc_cast(const float, (b ? (n - 1) : n)));
	}
}

MC_TARGET_FUNC void mc_mstdd1xnff(const int n, const float * x, const int b, double * mean, double * stddev)
{
	double sumsq, scale;

	*mean   = 0.0;
	*stddev = 0.0;
	if (n > 1) {
		mc_mssqr1xnff(n, x, (b > 1 ? 1 : 0), mean, &sumsq, &scale);
		*stddev = scale * mc_sqrt(sumsq / mc_cast(const double, (b ? (n - 1) : n)));
	}
}

MC_TARGET_FUNC void mc_mstdd1xn(const int n, const double * x, const int b, double * mean, double * stddev)
{
	double sumsq, scale;

	*mean   = 0.0;
	*stddev = 0.0;
	if (n > 1) {
		mc_mssqr1xn(n, x, (b > 1 ? 1 : 0), mean, &sumsq, &scale);
		*stddev = scale * mc_sqrt(sumsq / mc_cast(const double, (b ? (n - 1) : n)));
	}
}

MC_TARGET_FUNC void mc_mstdd1xnl(const int n, const long double * x, const int b, long double * mean, long double * stddev)
{
	long double sumsq, scale;

	*mean   = 0.0L;
	*stddev = 0.0L;
	if (n > 1) {
		mc_mssqr1xnl(n, x, (b > 1 ? 1 : 0), mean, &sumsq, &scale);
		*stddev = scale * mc_sqrtl(sumsq / mc_cast(const long double, (b ? (n - 1) : n)));
	}
}

#endif /* !MC_MSTDD1XN_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_stdd1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/numa/mc_mssqr1xn.h>

#ifndef MC_STDD1XN_H
#define MC_STDD1XN_H

#pragma mark - mc_stdd1xn -

MC_TARGET_FUNC void mc_stdd1xnf(const int n, const float * x, const int b, float * stddev)
{
	float mean, sumsq, scale;

	*stddev = 0.0f;
	if (n > 1) {
		mc_mssqr1xnf(n, x, 0, &mean, &sumsq, &scale);
		*stddev = scale * mc_sqrtf(sumsq / mc_cast(const float, (b ? (n - 1) : n)));
	}
}

MC_TARGET_FUNC void mc_stdd1xnff(const int n, const float * x, const int b, double * stddev)
{
	double mean, sumsq, scale;

	*stddev = 0.0;
	if (n > 1) {
		mc_mssqr1xnff(n, x, 0, &mean, &sumsq, &scale);
		*stddev = scale * mc_sqrt(sumsq / mc_cast(const double, (b ? (n - 1) : n)));
	}
}

MC_TARGET_FUNC void mc_stdd1xn(const int n, const double * x, const int b, double * stddev)
{
	double mean, sumsq, scale;

	*stddev = 0.0;
	if (n > 1) {
		mc_mssqr1xn(n, x, 0, &mean, &sumsq, &scale);
		*stddev = scale * mc_sqrt(sumsq / mc_cast(const double, (b ? (n - 1) : n)));
	}
}

MC_TARGET_FUNC void mc_stdd1xnl(const int n, const long double * x, const int b, long double * stddev)
{
	long double mean, sumsq, scale;

	*stddev = 0.0L;
	if (n > 1) {
		mc_mssqr1xnl(n, x, 0, &mean, &sumsq, &scale);
		*stddev = scale * mc_sqrtl(sumsq / mc_cast(const long double, (b ? (n - 1) : n)));
	}
}

#endif /* !MC_STDD1XN_H */

/* EOF */
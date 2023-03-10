//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mstddmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/numa/mc_mssqrmx1.h>

#ifndef MC_MSTDDMX1_H
#define MC_MSTDDMX1_H

#pragma mark - mc_mstddmx1 -

MC_TARGET_FUNC void mc_mstddmx1f(const int m, const int n, const int j, const float * a, const int b, float * mean, float * stddev)
{
	float sumsq, scale;

	*mean   = 0.0f;
	*stddev = 0.0f;
	if (n > 1) {
		mc_mssqrmx1f(m, n, j, a, (b > 1 ? 1 : 0), mean, &sumsq, &scale);
		*stddev = scale * mc_sqrtf(sumsq / mc_cast(const float, (b ? (m - 1) : m)));
	}
}

MC_TARGET_FUNC void mc_mstddmx1ff(const int m, const int n, const int j, const float * a, const int b, double * mean, double * stddev)
{
	double sumsq, scale;

	*mean   = 0.0;
	*stddev = 0.0;
	if (n > 1) {
		mc_mssqrmx1ff(m, n, j, a, (b > 1 ? 1 : 0), mean, &sumsq, &scale);
		*stddev = scale * mc_sqrt(sumsq / mc_cast(const double, (b ? (m - 1) : m)));
	}
}

MC_TARGET_FUNC void mc_mstddmx1(const int m, const int n, const int j, const double * a, const int b, double * mean, double * stddev)
{
	double sumsq, scale;

	*mean   = 0.0;
	*stddev = 0.0;
	if (n > 1) {
		mc_mssqrmx1(m, n, j, a, (b > 1 ? 1 : 0), mean, &sumsq, &scale);
		*stddev = scale * mc_sqrt(sumsq / mc_cast(const double, (b ? (m - 1) : m)));
	}
}

MC_TARGET_FUNC void mc_mstddmx1l(const int m, const int n, const int j, const long double * a, const int b, long double * mean, long double * stddev)
{
	long double sumsq, scale;

	*mean   = 0.0L;
	*stddev = 0.0L;
	if (n > 1) {
		mc_mssqrmx1l(m, n, j, a, (b > 1 ? 1 : 0), mean, &sumsq, &scale);
		*stddev = scale * mc_sqrtl(sumsq / mc_cast(const long double, (b ? (m - 1) : m)));
	}
}

#endif /* !MC_MSTDDMX1_H */

/* EOF */
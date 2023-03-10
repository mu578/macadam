//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_stddmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/numa/mc_ssqrmx1.h>

#ifndef MC_STDDMX1_H
#define MC_STDDMX1_H

#pragma mark - mc_stddmx1 -

MC_TARGET_FUNC float mc_stddmx1f(const int m, const int n, const int j, const float * a, const int b)
{
	float stddev, sumsq, scale;

	stddev = 0.0f;
	if (n > 1) {
		mc_ssqrmx1f(m, n, j, a, &sumsq, &scale);
		stddev = scale * mc_sqrtf(sumsq / mc_cast(const float, (b ? (m - 1) : m)));
	}
	return stddev;
}

MC_TARGET_FUNC double mc_stddmx1ff(const int m, const int n, const int j, const float * a, const int b)
{
	double stddev, sumsq, scale;

	stddev = 0.0;
	if (n > 1) {
		mc_ssqrmx1ff(m, n, j, a, &sumsq, &scale);
		stddev = scale * mc_sqrt(sumsq / mc_cast(const double, (b ? (m - 1) : m)));
	}
	return stddev;
}

MC_TARGET_FUNC double mc_stddmx1(const int m, const int n, const int j, const double * a, const int b)
{
	double stddev, sumsq, scale;

	stddev = 0.0;
	if (n > 1) {
		mc_ssqrmx1(m, n, j, a, &sumsq, &scale);
		stddev = scale * mc_sqrt(sumsq / mc_cast(const double, (b ? (m - 1) : m)));
	}
	return stddev;
}

MC_TARGET_FUNC long double mc_stddmx1l(const int m, const int n, const int j, const long double * a, const int b)
{
	long double stddev, sumsq, scale;

	stddev = 0.0L;
	if (n > 1) {
		mc_ssqrmx1l(m, n, j, a, &sumsq, &scale);
		stddev = scale * mc_sqrtl(sumsq / mc_cast(const long double, (b ? (m - 1) : m)));
	}
	return stddev;
}

#endif /* !MC_STDDMX1_H */

/* EOF */
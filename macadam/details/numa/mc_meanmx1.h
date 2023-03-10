//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_meanmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_summx1.h>

#ifndef MC_MEANMX1_H
#define MC_MEANMX1_H

#pragma mark - mc_meanmx1 -

MC_TARGET_FUNC float mc_meanmx1f(const int m, const int n, const int j, const float * a, const int b, const int f)
{
	float s = 0.0f;
	if (n > 0) {
		s = mc_summx1f(m, n, j, a, f);
		s = s / mc_cast(const float, (b ? m - 1 : m));
	}
	return s;
}

MC_TARGET_FUNC double mc_meanmx1ff(const int m, const int n, const int j, const float * a, const int b, const int f)
{
	double s = 0.0;
	if (n > 0) {
		s = mc_summx1ff(m, n, j, a, f);
		s = s / mc_cast(const double, (b ? m - 1 : m));
	}
	return s;
}

MC_TARGET_FUNC double mc_meanmx1(const int m, const int n, const int j, const double * a, const int b, const int f)
{
	double s = 0.0;
	if (n > 0) {
		s = mc_summx1(m, n, j, a, f);
		s = s / mc_cast(const double, (b ? m - 1 : m));
	}
	return s;
}

MC_TARGET_FUNC long double mc_meanmx1l(const int m, const int n, const int j, const long double * a, const int b, const int f)
{
	long double s = 0.0L;
	if (n > 0) {
		s = mc_summx1l(m, n, j, a, f);
		s = s / mc_cast(const long double, (b ? m - 1 : m));
	}
	return s;
}

#endif /* !MC_MEANMX1_H */

/* EOF */
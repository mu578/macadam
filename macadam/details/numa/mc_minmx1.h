//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_minmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_minmaxmx1.h>

#ifndef MC_MINMX1_H
#define MC_MINMX1_H

#pragma mark - mc_minmx1 -

MC_TARGET_FUNC float mc_minmx1f(const int m, const int n, const int j, const float * a)
{
	float min;
	mc_minmaxmx1f(m, n, j, a, &min, MC_NULLPTR, MC_NULLPTR, MC_NULLPTR);
	return min;
}

MC_TARGET_FUNC double mc_minmx1ff(const int m, const int n, const int j, const float * a)
{
	double min;
	mc_minmaxmx1ff(m, n, j, a, &min, MC_NULLPTR, MC_NULLPTR, MC_NULLPTR);
	return min;
}

MC_TARGET_FUNC double mc_minmx1(const int m, const int n, const int j, const double * a)
{
	double min;
	mc_minmaxmx1(m, n, j, a, &min, MC_NULLPTR, MC_NULLPTR, MC_NULLPTR);
	return min;
}

MC_TARGET_FUNC long double mc_minmx1l(const int m, const int n, const int j, const long double * a)
{
	long double min;
	mc_minmaxmx1l(m, n, j, a, &min, MC_NULLPTR, MC_NULLPTR, MC_NULLPTR);
	return min;
}

#endif /* !MC_MAXMX1_H */

/* EOF */
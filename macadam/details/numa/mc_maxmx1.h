//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_maxmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_minmaxmx1.h>

#ifndef MC_MAXMX1_H
#define MC_MAXMX1_H

#pragma mark - mc_maxmx1 -

MC_TARGET_FUNC float mc_maxmx1f(const int m, const int n, const int j, const float * a)
{
	float max;
	mc_minmaxmx1f(m, n, j, a, MC_NULLPTR, &max, MC_NULLPTR, MC_NULLPTR);
	return max;
}

MC_TARGET_FUNC double mc_maxmx1ff(const int m, const int n, const int j, const float * a)
{
	double max;
	mc_minmaxmx1ff(m, n, j, a, MC_NULLPTR, &max, MC_NULLPTR, MC_NULLPTR);
	return max;
}

MC_TARGET_FUNC double mc_maxmx1(const int m, const int n, const int j, const double * a)
{
	double max;
	mc_minmaxmx1(m, n, j, a, MC_NULLPTR, &max, MC_NULLPTR, MC_NULLPTR);
	return max;
}

MC_TARGET_FUNC long double mc_maxmx1l(const int m, const int n, const int j, const long double * a)
{
	long double max;
	mc_minmaxmx1l(m, n, j, a, MC_NULLPTR, &max, MC_NULLPTR, MC_NULLPTR);
	return max;
}

#endif /* !MC_MAXMX1_H */

/* EOF */
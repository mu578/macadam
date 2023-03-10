//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_argmaxmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_minmaxmx1.h>

#ifndef MC_ARGMAXMX1_H
#define MC_ARGMAXMX1_H

#pragma mark - mc_argmaxmx1 -

MC_TARGET_FUNC int mc_argmaxmx1f(const int m, const int n, const int j, const float * a)
{
	int q;
	float max;
	mc_minmaxmx1f(m, n, j, a, MC_NULLPTR, &max, MC_NULLPTR, &q);
	return q;
}

MC_TARGET_FUNC int mc_argmaxmx1(const int m, const int n, const int j, const double * a)
{
	int q;
	double max;
	mc_minmaxmx1(m, n, j, a, MC_NULLPTR, &max, MC_NULLPTR, &q);
	return q;
}

MC_TARGET_FUNC int mc_argmaxmx1l(const int m, const int n, const int j, const long double * a)
{
	int q;
	long double max;
	mc_minmaxmx1l(m, n, j, a, MC_NULLPTR, &max, MC_NULLPTR, &q);
	return q;
}

#endif /* !MC_ARGMAXMX1_H */

/* EOF */
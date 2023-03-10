//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rescalemx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_minmaxmx1.h>

#ifndef MC_RESCALEMX1_H
#define MC_RESCALEMX1_H

#pragma mark - mc_rescalemx1 -

MC_TARGET_FUNC void mc_rescalemx1f(const int m, const int n, const int j, float * c, float a, float b)
{
	int i = 0;
	float min, max, range;

	mc_minmaxmx1f(m, n, j, c, &min, &max, MC_NULLPTR, MC_NULLPTR);
	range = max - min;
	range = range != 0.0f ? 1.0f / range : 1.0f;
	for (; i < m; i++) {
		c[(n * i) + j] = a + ((c[(n * i) + j] - min) * range) * (b - a);
	}
}

MC_TARGET_FUNC void mc_rescalemx1(const int m, const int n, const int j, double * c, double a, double b)
{
	int i = 0;
	double min, max, range;

	mc_minmaxmx1(m, n, j, c, &min, &max, MC_NULLPTR, MC_NULLPTR);
	range = max - min;
	range = range != 0.0 ? 1.0 / range : 1.0;
	for (; i < m; i++) {
		c[(n * i) + j] = a + ((c[(n * i) + j] - min) * range) * (b - a);
	}
}

MC_TARGET_FUNC void mc_rescalemx1l(const int m, const int n, const int j, long double * c, long double a, long double b)
{
	int i = 0;
	long double min, max, range;

	mc_minmaxmx1l(m, n, j, c, &min, &max, MC_NULLPTR, MC_NULLPTR);
	range = max - min;
	range = range != 0.0L ? 1.0L / range : 1.0L;
	for (; i < m; i++) {
		c[(n * i) + j] = a + ((c[(n * i) + j] - min) * range) * (b - a);
	}
}

#endif /* !MC_RESCALEMX1_H */

/* EOF */
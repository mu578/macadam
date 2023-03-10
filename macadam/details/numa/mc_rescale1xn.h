//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rescale1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_minmax1xn.h>

#ifndef MC_RESCALE1XN_H
#define MC_RESCALE1XN_H

#pragma mark - mc_rescale1xn -

MC_TARGET_FUNC void mc_rescale1xnf(const int n, float * x, float a, float b)
{
	int i = 0;
	float min, max, range;

	mc_minmax1xnf(n, x, &min, &max, MC_NULLPTR, MC_NULLPTR);
	range = max - min;
	range = range != 0.0f ? 1.0f / range : 1.0f;
	for (; i < n; i++) {
		x[i] = a + ((x[i] - min) * range) * (b - a);
	}
}

MC_TARGET_FUNC void mc_rescale1xn(const int n, double * x, double a, double b)
{
	int i = 0;
	double min, max, range;

	mc_minmax1xn(n, x, &min, &max, MC_NULLPTR, MC_NULLPTR);
	range = max - min;
	range = range != 0.0 ? 1.0 / range : 1.0;
	for (; i < n; i++) {
		x[i] = a + ((x[i] - min) * range) * (b - a);
	}
}

MC_TARGET_FUNC void mc_rescale1xnl(const int n, long double * x, long double a, long double b)
{
	int i = 0;
	long double min, max, range;

	mc_minmax1xnl(n, x, &min, &max, MC_NULLPTR, MC_NULLPTR);
	range = max - min;
	range = range != 0.0L ? 1.0L / range : 1.0L;
	for (; i < n; i++) {
		x[i] = a + ((x[i] - min) * range) * (b - a);
	}
}

#endif /* !MC_RESCALE1XN_H */

/* EOF */
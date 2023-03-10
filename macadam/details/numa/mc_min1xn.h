//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_min1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_minmax1xn.h>

#ifndef MC_MIN1XN_H
#define MC_MIN1XN_H

#pragma mark - mc_min1xn -

MC_TARGET_FUNC float mc_min1xnf(const int n, const float * x)
{
	float min;
	mc_minmax1xnf(n, x, &min, MC_NULLPTR, MC_NULLPTR, MC_NULLPTR);
	return min;
}

MC_TARGET_FUNC double mc_min1xnff(const int n, const float * x)
{
	double min;
	mc_minmax1xnff(n, x, &min, MC_NULLPTR, MC_NULLPTR, MC_NULLPTR);
	return min;
}

MC_TARGET_FUNC double mc_min1xn(const int n, const double * x)
{
	double min;
	mc_minmax1xn(n, x, &min, MC_NULLPTR, MC_NULLPTR, MC_NULLPTR);
	return min;
}

MC_TARGET_FUNC long double mc_min1xnl(const int n, const long double * x)
{
	long double min;
	mc_minmax1xnl(n, x, &min, MC_NULLPTR, MC_NULLPTR, MC_NULLPTR);
	return min;
}

#endif /* !MC_MIN1XN_H */

/* EOF */
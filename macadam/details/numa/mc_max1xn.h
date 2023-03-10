//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_max1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_minmax1xn.h>

#ifndef MC_MAX1XN_H
#define MC_MAX1XN_H

#pragma mark - mc_max1xn -

MC_TARGET_FUNC float mc_max1xnf(const int n, const float * x)
{
	float max;
	mc_minmax1xnf(n, x, MC_NULLPTR, &max, MC_NULLPTR, MC_NULLPTR);
	return max;
}

MC_TARGET_FUNC double mc_max1xnff(const int n, const float * x)
{
	double max;
	mc_minmax1xnff(n, x, MC_NULLPTR, &max, MC_NULLPTR, MC_NULLPTR);
	return max;
}

MC_TARGET_FUNC double mc_max1xn(const int n, const double * x)
{
	double max;
	mc_minmax1xn(n, x, MC_NULLPTR, &max, MC_NULLPTR, MC_NULLPTR);
	return max;
}

MC_TARGET_FUNC long double mc_max1xnl(const int n, const long double * x)
{
	long double max;
	mc_minmax1xnl(n, x, MC_NULLPTR, &max, MC_NULLPTR, MC_NULLPTR);
	return max;
}

#endif /* !MC_MAX1XN_H */

/* EOF */
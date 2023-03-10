//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_argmin1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_minmax1xn.h>

#ifndef MC_ARGMIN1XN_H
#define MC_ARGMIN1XN_H

#pragma mark - mc_argmin1xn -

MC_TARGET_FUNC int mc_argmin1xnf(const int n, const float * x)
{
	int p;
	float min;
	mc_minmax1xnf(n, x, &min, MC_NULLPTR, &p, MC_NULLPTR);
	return p;
}

MC_TARGET_FUNC int mc_argmin1xn(const int n, const double * x)
{
	int p;
	double min;
	mc_minmax1xn(n, x, &min, MC_NULLPTR, &p, MC_NULLPTR);
	return p;
}

MC_TARGET_FUNC long double mc_argmin1xnl(const int n, const long double * x)
{
	int p;
	long double min;
	mc_minmax1xnl(n, x, &min, MC_NULLPTR, &p, MC_NULLPTR);
	return p;
}

#endif /* !MC_MIN1XN_H */

/* EOF */
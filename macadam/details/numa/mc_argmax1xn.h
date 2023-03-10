//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_argmax1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_minmax1xn.h>

#ifndef MC_ARGMAX1XN_H
#define MC_ARGMAX1XN_H

#pragma mark - mc_argmax1xn -

MC_TARGET_FUNC int mc_argmax1xnf(const int n, const float * x)
{
	int q;
	float max;
	mc_minmax1xnf(n, x, MC_NULLPTR, &max, MC_NULLPTR, &q);
	return q;
}

MC_TARGET_FUNC int mc_argmax1xn(const int n, const double * x)
{
	int q;
	double max;
	mc_minmax1xn(n, x, MC_NULLPTR, &max, MC_NULLPTR, &q);
	return q;
}

MC_TARGET_FUNC int mc_argmax1xnl(const int n, const long double * x)
{
	int q;
	long double max;
	mc_minmax1xnl(n, x, MC_NULLPTR, &max, MC_NULLPTR, &q);
	return q;
}

#endif /* !MC_ARGMAX1XN_H */

/* EOF */
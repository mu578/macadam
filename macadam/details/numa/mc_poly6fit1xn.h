//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_poly6fit1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_polyfit1xn.h>

#ifndef MC_POLY6FIT1XN_H
#define MC_POLY6FIT1XN_H

#pragma mark - mc_poly6fit1xn -

MC_TARGET_FUNC int mc_poly6fit1xnf(const int n, const float * x, const float * y, float c[7])
{
	float w[48];
	return mc_polyfit1xnf(n, 6, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly6fit1xnff(const int n, const float * x, const float * y, double c[7])
{
	double w[48];
	return mc_polyfit1xnff(n, 6, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly6fit1xn(const int n, const double * x, const double * y, double c[7])
{
	double w[48];
	return mc_polyfit1xn(n, 6, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly6fit1xnl(const int n, const long double * x, const long double * y, long double c[7])
{
	long double w[48];
	return mc_polyfit1xnl(n, 6, x, y, w, c);
}

#endif /* !MC_POLY6FIT1XN_H */

/* EOF */
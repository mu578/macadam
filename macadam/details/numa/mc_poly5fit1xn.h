//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_poly5fit1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_polyfit1xn.h>

#ifndef MC_POLY5FIT1XN_H
#define MC_POLY5FIT1XN_H

#pragma mark - mc_poly5fit1xn -

MC_TARGET_FUNC int mc_poly5fit1xnf(const int n, const float * x, const float * y, float c[6])
{
	float w[32];
	return mc_polyfit1xnf(n, 5, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly5fit1xnff(const int n, const float * x, const float * y, double c[6])
{
	double w[32];
	return mc_polyfit1xnff(n, 5, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly5fit1xn(const int n, const double * x, const double * y, double c[6])
{
	double w[32];
	return mc_polyfit1xn(n, 5, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly5fit1xnl(const int n, const long double * x, const long double * y, long double c[6])
{
	long double w[32];
	return mc_polyfit1xnl(n, 5, x, y, w, c);
}

#endif /* !MC_POLY5FIT1XN_H */

/* EOF */
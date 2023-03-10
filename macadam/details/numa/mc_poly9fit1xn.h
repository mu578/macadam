//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_poly9fit1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_polyfit1xn.h>

#ifndef MC_POLY9FIT1XN_H
#define MC_POLY9FIT1XN_H

#pragma mark - mc_poly9fit1xn -

MC_TARGET_FUNC int mc_poly9fit1xnf(const int n, const float * x, const float * y, float c[10])
{
	float w[64];
	return mc_polyfit1xnf(n, 9, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly9fit1xnff(const int n, const float * x, const float * y, double c[10])
{
	double w[64];
	return mc_polyfit1xnff(n, 9, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly9fit1xn(const int n, const double * x, const double * y, double c[10])
{
	double w[64];
	return mc_polyfit1xn(n, 9, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly9fit1xnl(const int n, const long double * x, const long double * y, long double c[10])
{
	long double w[64];
	return mc_polyfit1xnl(n, 9, x, y, w, c);
}

#endif /* !MC_POLY9FIT1XN_H */

/* EOF */
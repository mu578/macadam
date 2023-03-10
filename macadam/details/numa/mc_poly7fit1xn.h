//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_poly7fit1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_polyfit1xn.h>

#ifndef MC_POLY7FIT1XN_H
#define MC_POLY7FIT1XN_H

#pragma mark - mc_poly7fit1xn -

MC_TARGET_FUNC int mc_poly7fit1xnf(const int n, const float * x, const float * y, float c[8])
{
	float w[48];
	return mc_polyfit1xnf(n, 7, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly7fit1xnff(const int n, const float * x, const float * y, double c[8])
{
	double w[48];
	return mc_polyfit1xnff(n, 7, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly7fit1xn(const int n, const double * x, const double * y, double c[8])
{
	double w[48];
	return mc_polyfit1xn(n, 7, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly7fit1xnl(const int n, const long double * x, const long double * y, long double c[8])
{
	long double w[48];
	return mc_polyfit1xnl(n, 7, x, y, w, c);
}

#endif /* !MC_POLY7FIT1XN_H */

/* EOF */
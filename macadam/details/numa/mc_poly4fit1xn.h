//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_poly4fit1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_polyfit1xn.h>

#ifndef MC_POLY4FIT1XN_H
#define MC_POLY4FIT1XN_H

#pragma mark - mc_poly4fit1xn -

MC_TARGET_FUNC int mc_poly4fit1xnf(const int n, const float * x, const float * y, float c[5])
{
	float w[24];
	return mc_polyfit1xnf(n, 4, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly4fit1xnff(const int n, const float * x, const float * y, double c[5])
{
	double w[24];
	return mc_polyfit1xnff(n, 4, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly4fit1xn(const int n, const double * x, const double * y, double c[5])
{
	double w[24];
	return mc_polyfit1xn(n, 4, x, y, w, c);
}

MC_TARGET_FUNC int mc_poly4fit1xnl(const int n, const long double * x, const long double * y, long double c[5])
{
	long double w[24];
	return mc_polyfit1xnl(n, 4, x, y, w, c);
}

#endif /* !MC_POLY4FIT1XN_H */

/* EOF */
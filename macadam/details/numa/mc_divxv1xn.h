//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_divxv1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_DIVXV1XN_H
#define MC_DIVXV1XN_H

#pragma mark - mc_divxv1xn -

MC_TARGET_FUNC void mc_divxv1xnf(const int n, float * u, const float * x, float v)
{
//!# Ui=Xi/v
	const float w = 1.0f / v;

	int i = 0;
	for (; i < n; i++) {
		u[i] = x[i] * w;
	}
}

MC_TARGET_FUNC void mc_divxv1xnff(const int n, double * u, const float * x, float v)
{
//!# Ui=Xi/v
	const float w = 1.0f / v;

	int i = 0;
	for (; i < n; i++) {
		u[i] = mc_cast(double, x[i]) * mc_cast(double, w);
	}
}

MC_TARGET_FUNC void mc_divxv1xn(const int n, double * u, const double * x, double v)
{
//!# Ui=Xi/v
	const double w = 1.0 / v;

	int i = 0;
	for (; i < n; i++) {
		u[i] = x[i] * w;
	}
}

MC_TARGET_FUNC void mc_divxv1xnl(const int n, long double * u, const long double * x, long double v)
{
//!# Ui=Xi/v
	const long double w = 1.0L / v;

	int i = 0;
	for (; i < n; i++) {
		u[i] = x[i] * w;
	}
}

#endif /* !MC_DIVXV1XN_H */

/* EOF */
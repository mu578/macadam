//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mulxv1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_MULXV1XN_H
#define MC_MULXV1XN_H

#pragma mark - mc_mulxv1xn -

MC_TARGET_FUNC void mc_mulxv1xnf(const int n, float * u, const float * x, float v)
{
//!# Ui=Xi*v
	int i = 0;
	for (; i < n; i++) {
		u[i] = x[i] * v;
	}
}

MC_TARGET_FUNC void mc_mulxv1xnff(const int n, double * u, const float * x, float v)
{
//!# Ui=Xi*v
	int i = 0;
	for (; i < n; i++) {
		u[i] = mc_cast(double, x[i]) * mc_cast(double, v);
	}
}

MC_TARGET_FUNC void mc_mulxv1xn(const int n, double * u, const double * x, double v)
{
//!# Ui=Xi*v
	int i = 0;
	for (; i < n; i++) {
		u[i] = x[i] * v;
	}
}

MC_TARGET_FUNC void mc_mulxv1xnl(const int n, long double * u, const long double * x, long double v)
{
//!# Ui=Xi*v
	int i = 0;
	for (; i < n; i++) {
		u[i] = x[i] * v;
	}
}

#endif /* !MC_MULXV1XN_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_addxv1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_ADDXV1XN_H
#define MC_ADDXV1XN_H

#pragma mark - mc_addxv1xn -

MC_TARGET_FUNC void mc_addxv1xnf(const int n, float * u, const float * x, float v)
{
//!# Ui=Xi+v
	int i = 0;
	for (; i < n; i++) {
		u[i] = x[i] + v;
	}
}

MC_TARGET_FUNC void mc_addxv1xnff(const int n, double * u, const float * x, float v)
{
//!# Ui=Xi+v
	int i = 0;
	for (; i < n; i++) {
		u[i] = mc_cast(double, x[i]) * mc_cast(double, v);
	}
}

MC_TARGET_FUNC void mc_addxv1xn(const int n, double * u, const double * x, double v)
{
//!# Ui=Xi+v
	int i = 0;
	for (; i < n; i++) {
		u[i] = x[i] + v;
	}
}

MC_TARGET_FUNC void mc_addxv1xnl(const int n, long double * u, const long double * x, long double v)
{
//!# Ui=Xi+v
	int i = 0;
	for (; i < n; i++) {
		u[i] = x[i] + v;
	}
}

#endif /* !MC_ADDXV1XN_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_subxy1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_SUBXY1XN_H
#define MC_SUBXY1XN_H

#pragma mark - mc_subxy1xn -

MC_TARGET_FUNC void mc_subxy1xnf(const int n, float * u, const float * x, const float * y)
{
//!# Ui=Xi-Yi
	int i = 0;
	for (; i < n; i++) {
		u[i] = x[i] - y[i];
	}
}

MC_TARGET_FUNC void mc_subxy1xnff(const int n, double * u, const float * x, const float * y)
{
//!# Ui=Xi-Yi
	int i = 0;
	for (; i < n; i++) {
		u[i] = mc_cast(double, x[i]) * mc_cast(double, y[i]);
	}
}

MC_TARGET_FUNC void mc_subxy1xn(const int n, double * u, const double * x, const double * y)
{
//!# Ui=Xi-Yi
	int i = 0;
	for (; i < n; i++) {
		u[i] = x[i] - y[i];
	}
}

MC_TARGET_FUNC void mc_subxy1xnl(const int n, long double * u, const long double * x, const long double * y)
{
//!# Ui=Xi-Yi
	int i = 0;
	for (; i < n; i++) {
		u[i] = x[i] - y[i];
	}
}

#endif /* !MC_SUBXY1XN_H */

/* EOF */
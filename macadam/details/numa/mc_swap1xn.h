//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_swap1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcswap.h>

#ifndef MC_SWAP1XN_H
#define MC_SWAP1XN_H

#pragma mark - mc_swap1xn -

MC_TARGET_FUNC void mc_swap1xnf(const int n, float * MC_TARGET_RESTRICT x, float * MC_TARGET_RESTRICT y)
{
	int i = 0;
	float w;
	if (n > 0) {
		for (; i < n; i++) {
			mcswap_var(w, x[i], y[i]);
		}
	}
}

MC_TARGET_FUNC void mc_swap1xnfd(const int n, float * MC_TARGET_RESTRICT x, double * MC_TARGET_RESTRICT y)
{
	int i = 0;
	double w, xi, yi;
	if (n > 0) {
		for (; i < n; i++) {
			xi   = mc_cast(double, x[i]);
			yi   = y[i];
			mcswap_var(w, xi, yi);
			x[i] = mc_cast(float, yi);
			y[i] = xi;
		}
	}
}

MC_TARGET_FUNC void mc_swap1xndf(const int n, double * MC_TARGET_RESTRICT x, float * MC_TARGET_RESTRICT y)
{
	int i = 0;
	double w, xi, yi;
	if (n > 0) {
		for (; i < n; i++) {
			xi   = x[i];
			yi   = mc_cast(double, y[i]);
			mcswap_var(w, xi, yi);
			x[i] = yi;
			y[i] = mc_cast(float, xi);
		}
	}
}

MC_TARGET_FUNC void mc_swap1xn(const int n, double * MC_TARGET_RESTRICT x, double * MC_TARGET_RESTRICT y)
{
	int i = 0;
	double w;
	if (n > 0) {
		for (; i < n; i++) {
			mcswap_var(w, x[i], y[i]);
		}
	}
}

MC_TARGET_FUNC void mc_swap1xnl(const int n, long double * MC_TARGET_RESTRICT x, long double * MC_TARGET_RESTRICT y)
{
	int i = 0;
	long double w;
	if (n > 0) {
		for (; i < n; i++) {
			mcswap_var(w, x[i], y[i]);
		}
	}
}

#endif /* !MC_SWAP1XN_H */

/* EOF */
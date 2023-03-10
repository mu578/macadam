//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fliprg1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcswap.h>

#ifndef MC_FLIPRG1XN_H
#define MC_FLIPRG1XN_H

#pragma mark - mc_fliprg1xn -

MC_TARGET_FUNC void mc_fliprg1xnf(const int n, int p, int q, float * y, const float * x)
{
//!# Requires y[n] and x[n] where 1 < n. Y and X may be the same.
//!# Reversing order of elements in vector X and storing the result
//!# into vector Y within p and q range.
	int i = p, j = q;
	float w;
	if (p + q < n) {
		if (x != y) {
			for (; i < j; i++, j--) {
				y[i] = x[j];
			}
		} else {
			for (; i < j; i++, j--) {
				mcswap_var(w, y[i], y[j]);
			}
		}
	}
}

MC_TARGET_FUNC void mc_fliprg1xnff(const int n, int p, int q, double * y, const float * x)
{
//!# Requires y[n] and x[n] where 1 < n. Reversing order of elements in
//!# vector X and storing the result into vector Y within p and q range.
	int i = p, j = q;
	if (p + q < n) {
		for (; i < j; i++, j--) {
			y[i] = mc_cast(double, x[j]);
		}
	}
}

MC_TARGET_FUNC void mc_fliprg1xn(const int n, int p, int q, double * y, const double * x)
{
//!# Requires y[n] and x[n] where 1 < n. Y and X may be the same.
//!# Reversing order of elements in vector X and storing the result
//!# into vector Y within p and q range.
	int i = p, j = q;
	double w;
	if (p + q < n) {
		if (x != y) {
			for (; i < j; i++, j--) {
				y[i] = x[j];
			}
		} else {
			for (; i < j; i++, j--) {
				mcswap_var(w, y[i], y[j]);
			}
		}
	}
}

MC_TARGET_FUNC void mc_fliprg1xnl(const int n, int p, int q, long double * y, const long double * x)
{
//!# Requires y[n] and x[n] where 1 < n. Y and X may be the same.
//!# Reversing order of elements in vector X and storing the result
//!# into vector Y within p and q range.
	int i = p, j = q;
	long double w;
	if (p + q < n) {
		if (x != y) {
			for (; i < j; i++, j--) {
				y[i] = x[j];
			}
		} else {
			for (; i < j; i++, j--) {
				mcswap_var(w, y[i], y[j]);
			}
		}
	}
}

#endif /* !MC_FLIPRG1XN_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rotl1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_gcd.h>

#ifndef MC_ROTL1XN_H
#define MC_ROTL1XN_H

#pragma mark - mc_rotl1xn -

MC_TARGET_FUNC void mc_rotl1xnf(const int n, const int k, float * y, const float * x)
{
//!# Requires y[n] and x[n] where 1 < n. Y and X may be the same.
//!# Left rotate elements of vector X by k times and store the
//!# result into vector Y.
	int i = 0, j, l; 
	float w;
	if (x != y) {
		mc_copy1xnf(n, y , x);
	}
	for (; i < mc_igcd(k, n); i++) {
		j = i;
		w = y[j];
		while (!(l == i)) { 
			l    = ((j + k) >= n) ? ((j + k) - n) : (j + k);
			y[j] = y[l]; 
			j    = l; 
		}
		y[j] = w; 
	}
}

MC_TARGET_FUNC void mc_rotl1xnff(const int n, const int k, double * y, const float * x)
{
//!# Requires y[n] and x[n] where 1 < n. Left rotate elements of
//!# vector X by k times and store the result into vector Y.
	int i = 0, j, l; 
	double w;
	mc_copy1xnff(n, y , x);
	for (; i < mc_igcd(k, n); i++) {
		j = i;
		w = y[j];
		while (!(l == i)) { 
			l    = ((j + k) >= n) ? ((j + k) - n) : (j + k);
			y[j] = y[l]; 
			j    = l; 
		}
		y[j] = w; 
	}
}

MC_TARGET_FUNC void mc_rotl1xn(const int n, const int k, double * y, const double * x)
{
//!# Requires y[n] and x[n] where 1 < n. Y and X may be the same.
//!# Left rotate elements of vector X by k times and store the
//!# result into vector Y.
	int i = 0, j, l; 
	double w;
	if (x != y) {
		mc_copy1xn(n, y , x);
	}
	for (; i < mc_igcd(k, n); i++) {
		j = i;
		w = y[j];
		while (!(l == i)) { 
			l    = ((j + k) >= n) ? ((j + k) - n) : (j + k);
			y[j] = y[l]; 
			j    = l; 
		}
		y[j] = w; 
	}
}

MC_TARGET_FUNC void mc_rotl1xnl(const int n, const int k, long double * y, const long double * x)
{
//!# Requires y[n] and x[n] where 1 < n. Y and X may be the same.
//!# Left rotate elements of vector X by k times and store the
//!# result into vector Y.
	int i = 0, j, l; 
	long double w;
	if (x != y) {
		mc_copy1xnl(n, y , x);
	}
	for (; i < mc_igcd(k, n); i++) {
		j = i;
		w = y[j];
		while (!(l == i)) { 
			l    = ((j + k) >= n) ? ((j + k) - n) : (j + k);
			y[j] = y[l]; 
			j    = l; 
		}
		y[j] = w; 
	}
}

#endif /* !MC_ROTLR1XN_H */

/* EOF */
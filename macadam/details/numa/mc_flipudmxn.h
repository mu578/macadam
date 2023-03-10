//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_flipudmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_flipmx1.h>

#ifndef MC_FLIPUDMXN_H
#define MC_FLIPUDMXN_H

#pragma mark - mc_flipudmxn -

MC_TARGET_FUNC void mc_flipudmxnf(const int m, const int n, float * c, const float * a)
{
//!# Requires c[m x n] and a[m x n]. C and A may be the same. Reversing
//!# order of elements in each column of A and storing the result into C
//!# i.e returns A with its rows flipped in the up-down direction.
	int j = 0;
	for (; j < n; j++) {
		mc_flipmx1f(m, n, j, c, a);
	}
}

MC_TARGET_FUNC void mc_flipudmxnff(const int m, const int n, double * c, const float * a)
{
//!# Requires c[m x n] and a[m x n]. Reversing order of elements in each
//!# column of A and storing the result into C i.e returns A with its rows
//!# flipped in the up-down direction.
	int j = 0;
	for (; j < n; j++) {
		mc_flipmx1ff(m, n, j, c, a);
	}
}

MC_TARGET_FUNC void mc_flipudmxn(const int m, const int n, double * c, const double * a)
{
//!# Requires c[m x n] and a[m x n]. C and A may be the same. Reversing
//!# order of elements in each column of A and storing the result into C
//!# i.e returns A with its rows flipped in the up-down direction.
	int j = 0;
	for (; j < n; j++) {
		mc_flipmx1(m, n, j, c, a);
	}
}

MC_TARGET_FUNC void mc_flipudmxnl(const int m, const int n, long double * c, const long double * a)
{
//!# Requires c[m x n] and a[m x n]. C and A may be the same. Reversing
//!# order of elements in each column of A and storing the result into C
//!# i.e returns A with its rows flipped in the up-down direction.
	int j = 0;
	for (; j < n; j++) {
		mc_flipmx1l(m, n, j, c, a);
	}
}

#endif /* !MC_FLIPUDMXN_H */

/* EOF */
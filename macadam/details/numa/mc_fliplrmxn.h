//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_fliplrmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_flip1xn.h>

#ifndef MC_FLIPLRMXN_H
#define MC_FLIPLRMXN_H

#pragma mark - mc_fliplrmxn -

MC_TARGET_FUNC void mc_fliplrmxnf(const int m, const int n, float * c, const float * a)
{
//!# Requires c[m x n] and a[m x n]. C and A may be the same. Reversing
//!# order of elements in each row of A and storing the result into C
//!# i.e returns A with its columns flipped in the left-right direction.
	int i = 0;
	if (a != c) {
		for (; i < m; i++) {
			const float * x = a + (n * i);
			float * y       = c + (n * i);
			mc_flip1xnf(n, y, x);
		}
	} else {
		for (; i < m; i++) {
			float * y = c + (n * i);
			mc_flip1xnf(n, y, y);
		}
	}
}

MC_TARGET_FUNC void mc_fliplrmxnff(const int m, const int n, double * c, const float * a)
{
//!# Requires c[m x n] and a[m x n]. Reversing order of elements in each
//!# row of A and storing the result into C i.e returns A with its columns
//!# flipped in the left-right direction.
	int i = 0;
	for (; i < m; i++) {
		const float  * x = a + (n * i);
		double * y       = c + (n * i);
		mc_flip1xnff(n, y, x);
	}
}

MC_TARGET_FUNC void mc_fliplrmxn(const int m, const int n, double * c, const double * a)
{
//!# Requires c[m x n] and a[m x n]. C and A may be the same. Reversing
//!# order of elements in each row of A and storing the result into C
//!# i.e returns A with its columns flipped in the left-right direction.
	int i = 0;
	if (a != c) {
		for (; i < m; i++) {
			const double * x = a + (n * i);
			double * y       = c + (n * i);
			mc_flip1xn(n, y, x);
		}
	} else {
		for (; i < m; i++) {
			double * y = c + (n * i);
			mc_flip1xn(n, y, y);
		}
	}
}

MC_TARGET_FUNC void mc_fliplrmxnl(const int m, const int n, long double * c, const long double * a)
{
//!# Requires c[m x n] and a[m x n]. C and A may be the same. Reversing
//!# order of elements in each row of A and storing the result into C
//!# i.e returns A with its columns flipped in the left-right direction.
	int i = 0;
	if (a != c) {
		for (; i < m; i++) {
			const long double * x = a + (n * i);
			long double * y       = c + (n * i);
			mc_flip1xnl(n, y, x);
		}
	} else {
		for (; i < m; i++) {
			long double * y = c + (n * i);
			mc_flip1xnl(n, y, y);
		}
	}
}

#endif /* !MC_FLIPLRMXN_H */

/* EOF */
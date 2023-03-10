//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_trsinxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_TRSISNXN_H
#define MC_TRSISNXN_H

#pragma mark - mc_trsinxn -

MC_TARGET_FUNC void mc_trsinxnf(const int n, float * a)
{
//!# Returning transpose of A (in-place).
	int i = 0, j;
	float w;
	for (; i < (n - 2); i++) {
		for (j = (i + 1); j < (n - 1); j++) {
			w                = a[((i) * n) + j];
			a[((i) * n) + j] = a[((j) * n) + i];
			a[((j) * n) + i] = w;
		}
	}
}

MC_TARGET_FUNC void mc_trsinxn(const int n, double * a)
{
//!# Returning transpose of A (in-place).
	int i = 0, j;
	double w;
	for (; i < (n - 2); i++) {
		for (j = (i + 1); j < (n - 1); j++) {
			w                = a[((i) * n) + j];
			a[((i) * n) + j] = a[((j) * n) + i];
			a[((j) * n) + i] = w;
		}
	}
}

MC_TARGET_FUNC void mc_trsinxnl(const int n, long double * a)
{
//!# Returning transpose of A (in-place).
	int i = 0, j;
	long double w;
	for (; i < (n - 2); i++) {
		for (j = (i + 1); j < (n - 1); j++) {
			w                = a[((i) * n) + j];
			a[((i) * n) + j] = a[((j) * n) + i];
			a[((j) * n) + i] = w;
		}
	}
}

#endif /* !MC_TRSISNXN_H */

/* EOF */
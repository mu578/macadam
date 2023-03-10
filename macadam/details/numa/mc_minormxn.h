//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_minormxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_trsimxn.h>

#ifndef MC_MINORMXN_H
#define MC_MINORMXN_H

#pragma mark - mc_minormxn -

MC_TARGET_FUNC void mc_minormxnf(const int m, const int n, const int p, float * MC_TARGET_RESTRICT b, const float * a)
{
	int i = 0, j;
	for (; i < p; i++) {
		b[(n * i) + i] = 1.0f;
	}
	for (i = p; i < m; i++) {
		for (j = p; j < n; j++) {
			b[(n * i) + j] = a[(n * i) + j];
		}
	}
}

MC_TARGET_FUNC void mc_minormxnff(const int m, const int n, const int p, double * b, const float * a)
{
	int i = 0, j;
	for (; i < p; i++) {
		b[(n * i) + i] = 1.0;
	}
	for (i = p; i < m; i++) {
		for (j = p; j < n; j++) {
			b[(n * i) + j] = mc_cast(double, a[(n * i) + j]);
		}
	}
}

MC_TARGET_FUNC void mc_minormxn(const int m, const int n, const int p, double * b, const double * a)
{
	int i = 0, j;
	for (; i < p; i++) {
		b[(n * i) + i] = 1.0;
	}
	for (i = p; i < m; i++) {
		for (j = p; j < n; j++) {
			b[(n * i) + j] = a[(n * i) + j];
		}
	}
}

MC_TARGET_FUNC void mc_minormxnl(const int m, const int n, const int p, long double * MC_TARGET_RESTRICT b, const long double * a)
{
	int i = 0, j;
	for (; i < p; i++) {
		b[(n * i) + i] = 1.0;
	}
	for (i = p; i < m; i++) {
		for (j = p; j < n; j++) {
			b[(n * i) + j] = a[(n * i) + j];
		}
	}
}

#endif /* !MC_MINORMXN_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_nnzmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_NNZMX1_H
#define MC_NNZMX1_H

#pragma mark - mc_nnzmx1 -

MC_TARGET_FUNC int mc_nnzmx1f(const int m, const int n, const int j, const float * a)
{
	int i = 0, c = 0;
	for (; i < m; i++) {
		if (a[(n * i) + j] != 0.0f) { ++c; }
	}
	return c;
}

MC_TARGET_FUNC int mc_nnzmx1(const int m, const int n, const int j, const double * a)
{
	int i = 0, c = 0;
	for (; i < m; i++) {
		if (a[(n * i) + j] != 0.0) { ++c; }
	}
	return c;
}

MC_TARGET_FUNC int mc_nnzmx1l(const int m, const int n, const int j, const long double * a)
{
	int i = 0, c = 0;
	for (; i < m; i++) {
		if (a[(n * i) + j] != 0.0L) { ++c; }
	}
	return c;
}

#endif /* !MC_NNZMX1_H */

/* EOF */
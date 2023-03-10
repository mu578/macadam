//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_unitmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_l2normmx1.h>

#ifndef MC_UNITMX1_H
#define MC_UNITMX1_H

#pragma mark - mc_unitmx1 -

MC_TARGET_FUNC void mc_unitmx1f(const int m, const int n, const int j, float * a)
{
	const float norm = mc_l2normmx1f(m, n, j, a);
	if (norm != 0.0f) {
		const float scale = 1.0f / norm;
		int i             = 0;
		for (; i < m; i++) {
			a[(n * i) + j] = a[(n * i) + j] * scale;
		}
	} else {
		a[(n * j) + j] = 1.0f;
	}
}

MC_TARGET_FUNC void mc_unitmx1(const int m, const int n, const int j, double * a)
{
	const double norm = mc_l2normmx1(m, n, j, a);
	if (norm != 0.0) {
		const double scale = 1.0 / norm;
		int i              = 0;
		for (; i < m; i++) {
			a[(n * i) + j] = a[(n * i) + j] * scale;
		}
	} else {
		a[(n * j) + j] = 1.0;
	}
}

MC_TARGET_FUNC void mc_unitmx1l(const int m, const int n, const int j, long double * a)
{
	const long double norm = mc_l2normmx1l(m, n, j, a);
	if (norm != 0.0L) {
		const long double scale = 1.0L / norm;
		int i                   = 0;
		for (; i < m; i++) {
			a[(n * i) + j] = a[(n * i) + j] * scale;
		}
	} else {
		a[(n * j) + j] = 1.0L;
	}
}

#endif /* !MC_UNITMX1_H */

/* EOF */
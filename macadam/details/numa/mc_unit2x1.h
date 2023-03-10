//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_unit2x1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_l2norm2x1.h>

#ifndef MC_UNIT2X1_H
#define MC_UNIT2X1_H

#pragma mark - mc_unit2x1 -

MC_TARGET_FUNC void mc_unit2x1f(const int n, const int j, float * a)
{
	const float norm = mc_l2norm2x1f(n, j, a);
	if (norm != 0.0f) {
		const float scale = 1.0f / norm;
		a[j]              = a[j] * scale;
		a[n + j]          = a[n + j] * scale;
	} else {
		a[(n * j) + j] = 1.0f;
	}
}

MC_TARGET_FUNC void mc_unit2x1(const int n, const int j, double * a)
{
	const double norm = mc_l2norm2x1(n, j, a);
	if (norm != 0.0) {
		const double scale = 1.0 / norm;
		a[j]              = a[j] * scale;
		a[n + j]          = a[n + j] * scale;
	} else {
		a[(n * j) + j] = 1.0;
	}
}

MC_TARGET_FUNC void mc_unit2x1l(const int n, const int j, long double * a)
{
	const long double norm = mc_l2norm2x1l(n, j, a);
	if (norm != 0.0L) {
		const long double scale = 1.0L / norm;
		a[j]                    = a[j] * scale;
		a[n + j]                = a[n + j] * scale;
	} else {
		a[(n * j) + j] = 1.0L;
	}
}

#endif /* !MC_UNIT2X1_H */

/* EOF */
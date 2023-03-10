//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zeros3x1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_ZEROS3X1_H
#define MC_ZEROS3X1_H

#pragma mark - mc_zeros3x1 -

MC_TARGET_FUNC void mc_zeros3x1f(const int n, const int j, float * a)
{
	a[j]           = 0.0f;
	a[n + j]       = 0.0f;
	a[(n * 2) + j] = 0.0f;
}

MC_TARGET_FUNC void mc_zeros3x1(const int n, const int j, double * a)
{
	a[j]           = 0.0;
	a[n + j]       = 0.0;
	a[(n * 2) + j] = 0.0;
}

MC_TARGET_FUNC void mc_zeros3x1l(const int n, const int j, long double * a)
{
	a[j]           = 0.0L;
	a[n + j]       = 0.0L;
	a[(n * 2) + j] = 0.0L;
}

#endif /* !MC_ZEROS3X1_H */

/* EOF */
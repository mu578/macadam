//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_copy1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/details/mc_mem.h>

#ifndef MC_COPY1XN_H
#define MC_COPY1XN_H

#pragma mark - mc_copy1xn -

MC_TARGET_FUNC void mc_copy1xnf(const int n, float * y, const float * x)
{
	mc_base_memcpy_type(float, n, y, x);
}

MC_TARGET_FUNC void mc_copy1xnff(const int n, double * y, const float * x)
{
	int i = 0;
	for (; i < n; i++) {
		y[i] = mc_cast(double, x[i]);
	}
}

MC_TARGET_FUNC void mc_copy1xn(const int n, double * y, const double * x)
{
	mc_base_memcpy_type(double, n, y, x);
}

MC_TARGET_FUNC void mc_copy1xnl(const int n, long double * y, const long double * x)
{
	mc_base_memcpy_type(long double, n, y, x);
}

#endif /* !MC_COPY1XN_H */

/* EOF */
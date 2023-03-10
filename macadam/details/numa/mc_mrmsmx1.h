//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mrmsmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_mstddmx1.h>

#ifndef MC_MRMSMX1_H
#define MC_MRMSMX1_H

#pragma mark - mc_mrmsmx1 -

MC_TARGET_FUNC void mc_mrmsmx1f(const int m, const int n, const int j, const float * a, float * mean, float * rms)
{
	mc_mstddmx1f(m, n, j, a, 1, mean, rms);
}

MC_TARGET_FUNC void mc_mrmsmx1ff(const int m, const int n, const int j, const float * a, double * mean, double * rms)
{
	mc_mstddmx1ff(m, n, j, a, 1, mean, rms);
}

MC_TARGET_FUNC void mc_mrmsmx1(const int m, const int n, const int j, const double * a, double * mean, double * rms)
{
	mc_mstddmx1(m, n, j, a, 1, mean, rms);
}

MC_TARGET_FUNC void mc_mrmsmx1l(const int m, const int n, const int j, const long double * a, long double * mean, long double * rms)
{
	mc_mstddmx1l(m, n, j, a, 1, mean, rms);
}

#endif /* !MC_MRMSMX1_H */

/* EOF */
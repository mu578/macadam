//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rmsmx1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_stddmx1.h>

#ifndef MC_RMSMX1_H
#define MC_RMSMX1_H

#pragma mark - mc_rmsmx1 -

MC_TARGET_FUNC float mc_rmsmx1f(const int m, const int n, const int j, const float * a)
{
	return mc_stddmx1f(m, n, j, a, 0);
}

MC_TARGET_FUNC double mc_rmsmx1ff(const int m, const int n, const int j, const float * a)
{
	return mc_stddmx1ff(m, n, j, a, 0);
}

MC_TARGET_FUNC double mc_rmsmx1(const int m, const int n, const int j, const double * a)
{
	return mc_stddmx1(m, n, j, a, 0);
}

MC_TARGET_FUNC long double mc_rmsmx1l(const int m, const int n, const int j, const long double * a)
{
	return mc_stddmx1l(m, n, j, a, 0);
}

#endif /* !MC_RMSMX1_H */

/* EOF */
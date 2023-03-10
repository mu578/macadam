//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_asum1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_sum1xn.h>

#ifndef MC_ASUM1XN_H
#define MC_ASUM1XN_H

#pragma mark - mc_asum1xn -

MC_TARGET_FUNC float mc_asum1xnf(const int n, const float * x)
{
	return mc_sum1xnf(n, x, 4);
}

MC_TARGET_FUNC double mc_asum1xnff(const int n, const float * x)
{
	return mc_sum1xnff(n, x, 4);
}

MC_TARGET_FUNC double mc_asum1xn(const int n, const double * x)
{
	return mc_sum1xn(n, x, 4);
}

MC_TARGET_FUNC long double mc_asum1xnl(const int n, const long double * x)
{
	return mc_sum1xnl(n, x, 4);
}

#endif /* !MC_ASUM1XN_H */

/* EOF */
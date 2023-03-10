//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_l1norm1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_a2sum1xn.h>

#ifndef MC_L1NORM1XN_H
#define MC_L1NORM1XN_H

#pragma mark - mc_l1norm1xn -

MC_TARGET_FUNC float mc_l1norm1xnf(const int n, const float * x)
{
//!# Requires x[1 x n]. Returning the l1-norm of x.
	return mc_a2sum1xnf(n, x);
}

MC_TARGET_FUNC double mc_l1norm1xnff(const int n, const float * x)
{
//!# Requires x[1 x n]. Returning the l1-norm of x.
	return mc_a2sum1xnff(n, x);
}

MC_TARGET_FUNC double mc_l1norm1xn(const int n, const double * x)
{
//!# Requires x[1 x n]. Returning the l1-norm of x.
	return mc_a2sum1xn(n, x);
}

MC_TARGET_FUNC long double mc_l1norm1xnl(const int n, const long double * x)
{
//!# Requires x[1 x n]. Returning the l1-norm of x.
	return mc_a2sum1xnl(n, x);
}

#endif /* !MC_L1NORM1XN_H */

/* EOF */
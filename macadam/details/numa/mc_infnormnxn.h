//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_infnormnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_infnormmxn.h>

#ifndef MC_INFNORMNXN_H
#define MC_INFNORMNXN_H

#pragma mark - mc_infnormnxn -

MC_TARGET_FUNC float mc_infnormnxnf(const int n, const float * a, const int f)
{
//!# Requires a[n x n]. Returning the infinity norm of a.
//!# f=0: computing the maximum of the absolute row sums.
//!# f=1: computing the minimum of the absolute row sums.
	return mc_infnormmxnf(n, n, a, f);
}

MC_TARGET_FUNC double mc_infnormnxnff(const int n, const float * a, const int f)
{
//!# Requires a[n x n]. Returning the infinity norm of a.
//!# f=0: computing the maximum of the absolute row sums.
//!# f=1: computing the minimum of the absolute row sums.
	return mc_infnormmxnff(n, n, a, f);
}

MC_TARGET_FUNC double mc_infnormnxn(const int n, const double * a, const int f)
{
//!# Requires a[n x n]. Returning the infinity norm of a.
//!# f=0: computing the maximum of the absolute row sums.
//!# f=1: computing the minimum of the absolute row sums.
	return mc_infnormmxn(n, n, a, f);
}

MC_TARGET_FUNC long double mc_infnormnxnl(const int n, const long double * a, const int f)
{
//!# Requires a[n x n]. Returning the infinity norm of a.
//!# f=0: computing the maximum of the absolute row sums.
//!# f=1: computing the minimum of the absolute row sums.
	return mc_infnormmxnl(n, n, a, f);
}

#endif /* !MC_INFNORMNXN_H */

/* EOF */
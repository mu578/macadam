//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mulaxnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_mulaxmxn.h>

#ifndef MC_MULAXNXN_H
#define MC_MULAXNXN_H

#pragma mark - mc_mulaxnxn -

MC_TARGET_FUNC void mc_mulaxnxnf(const int n, float * MC_TARGET_RESTRICT b, const float * a, const float * x)
{
//!# Requires b[n x 1], a[n * n] and x[n x 1].
//!# b=a*x
	mc_mulaxmxnf(n, n, b, a, x);
}

MC_TARGET_FUNC void mc_mulaxnxnff(const int n, double * MC_TARGET_RESTRICT b, const float * a, const float * x)
{
//!# Requires b[n x 1], a[n * n] and x[n x 1].
//!# b=a*x
	mc_mulaxmxnff(n, n, b, a, x);
}

MC_TARGET_FUNC void mc_mulaxnxnfd(const int n, double * MC_TARGET_RESTRICT b, const float * a, const double * x)
{
//!# Requires b[n x 1], a[n * n] and x[n x 1].
//!# b=a*x
	mc_mulaxmxnfd(n, n, b, a, x);
}

MC_TARGET_FUNC void mc_mulaxnxndf(const int n, double * MC_TARGET_RESTRICT b, const double * a, const float * x)
{
//!# Requires b[n x 1], a[n * n] and x[n x 1].
//!# b=a*x
	mc_mulaxmxndf(n, n, b, a, x);
}

MC_TARGET_FUNC void mc_mulaxnxn(const int n, double * MC_TARGET_RESTRICT b, const double * a, const double * x)
{
//!# Requires b[n x 1], a[n * n] and x[n x 1].
//!# b=a*x
	mc_mulaxmxn(n, n, b, a, x);
}

MC_TARGET_FUNC void mc_mulaxnxnl(const int n, long double * MC_TARGET_RESTRICT b, const long double * a, const long double * x)
{
//!# Requires b[n x 1], a[n * n] and x[n x 1].
//!# b=a*x
	mc_mulaxmxnl(n, n, b, a, x);
}

#endif /* !MC_MULAXNXN_H */

/* EOF */
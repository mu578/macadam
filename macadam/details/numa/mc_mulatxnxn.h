//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mulatxnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_mulatxmxn.h>

#ifndef MC_MULATXNXN_H
#define MC_MULATXNXN_H

#pragma mark - mc_mulatxnxn -

MC_TARGET_FUNC void mc_mulatxnxnf(const int n, float * MC_TARGET_RESTRICT b, const float * a, const float * x)
{
//!# Requires b[n x 1], a[n * n] and x[n x 1].
//!# b=a'*x
	mc_mulatxmxnf(n, n, b, a, x);
}

MC_TARGET_FUNC void mc_mulatxnxnff(const int n, double * MC_TARGET_RESTRICT b, const float * a, const float * x)
{
//!# Requires b[n x 1], a[n * n] and x[n x 1].
//!# b=a'*x
	mc_mulatxmxnff(n, n, b, a, x);
}

MC_TARGET_FUNC void mc_mulatxnxnfd(const int n, double * MC_TARGET_RESTRICT b, const float * a, const double * x)
{
//!# Requires b[n x 1], a[n * n] and x[n x 1].
//!# b=a'*x
	mc_mulatxmxnfd(n, n, b, a, x);
}

MC_TARGET_FUNC void mc_mulatxnxndf(const int n, double * MC_TARGET_RESTRICT b, const double * a, const float * x)
{
//!# Requires b[n x 1], a[n * n] and x[n x 1].
//!# b=a'*x
	mc_mulatxmxndf(n, n, b, a, x);
}

MC_TARGET_FUNC void mc_mulatxnxn(const int n, double * MC_TARGET_RESTRICT b, const double * a, const double * x)
{
//!# Requires b[n x 1], a[n * n] and x[n x 1].
//!# b=a'*x
	mc_mulatxmxn(n, n, b, a, x);
}

MC_TARGET_FUNC void mc_mulatxnxnl(const int n, long double * MC_TARGET_RESTRICT b, const long double * a, const long double * x)
{
//!# Requires b[n x 1], a[n * n] and x[n x 1].
//!# b=a'*x
	mc_mulatxmxnl(n, n, b, a, x);
}

#endif /* !MC_MULATXNXN_H */

/* EOF */
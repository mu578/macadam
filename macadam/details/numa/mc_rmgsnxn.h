//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rmgsnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_orthrnxn.h>

#ifndef MC_RMGSNXN_H
#define MC_RMGSNXN_H

#pragma mark - mc_rmgsnxn -

MC_TARGET_FUNC int mc_rmgsnxnf(const int n, const float * a, float * q, float * MC_TARGET_RESTRICT r, float * MC_TARGET_RESTRICT w, int * pvi)
{
//!# Requires a[m x n], q[m x n] and r[n x n] where 1 < n.
//!# A and Q may be the same. Performing a QR decomposition using
//!# Recursive Modified Gram-Schmidt method with optional pivoting.
	return mc_orthrnxnf(n, a, MCLIMITS_TINYF, q, r , w, pvi);
}

MC_TARGET_FUNC int mc_rmgsnxnff(const int n, const float * a, double * q, double * MC_TARGET_RESTRICT r, double * MC_TARGET_RESTRICT w, int * pvi)
{
//!# Requires a[m x n], q[m x n] and r[n x n] where 1 < n.
//!# Performing a QR decomposition using Recursive Modified Gram-Schmidt
//!# method with optional pivoting.
	return mc_orthrnxnff(n, a, MCLIMITS_TINYF, q, r , w, pvi);
}

MC_TARGET_FUNC int mc_rmgsnxn(const int n, const double * a, double * q, double * MC_TARGET_RESTRICT r, double * MC_TARGET_RESTRICT w, int * pvi)
{
//!# Requires a[m x n], q[m x n] and r[n x n] where 1 < n.
//!# A and Q may be the same. Performing a QR decomposition using
//!# Recursive Modified Gram-Schmidt method with optional pivoting.
	return mc_orthrnxn(n, a, MCLIMITS_TINY, q, r , w, pvi);
}

MC_TARGET_FUNC int mc_rmgsnxnl(const int n, const long double * a, long double * q, long double * MC_TARGET_RESTRICT r, long double * MC_TARGET_RESTRICT w, int * pvi)
{
//!# Requires a[m x n], q[m x n] and r[n x n] where 1 < n.
//!# A and Q may be the same. Performing a QR decomposition using
//!# Recursive Modified Gram-Schmidt method with optional pivoting.
	return mc_orthrnxnl(n, a, MCLIMITS_TINYL, q, r , w, pvi);
}

#endif /* !MC_RMGSNXN_H */

/* EOF */
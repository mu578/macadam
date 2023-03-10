//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rmgsmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_orthrmxn.h>

#ifndef MC_RMGSMXN_H
#define MC_RMGSMXN_H

#pragma mark - mc_rmgsmxn -

MC_TARGET_FUNC int mc_rmgsmxnf(const int m, const int n, const float * a, float * q, float * MC_TARGET_RESTRICT r, float * MC_TARGET_RESTRICT w, int * pvi)
{
//!# Requires a[m x n], q[m x n] and r[n x n] where 1 < n <= m.
//!# A and Q may be the same. Performing a QR decomposition using
//!# Recursive Modified Gram-Schmidt method with optional pivoting.
	return mc_orthrmxnf(m, n, a, MCLIMITS_TINYF, q, r, w, pvi);
}

MC_TARGET_FUNC int mc_rmgsmxnff(const int m, const int n, const float * a, double * q, double * MC_TARGET_RESTRICT r, double * MC_TARGET_RESTRICT w, int * pvi)
{
//!# Requires a[m x n], q[m x n] and r[n x n] where 1 < n <= m.
//!# Performing a QR decomposition using Recursive Modified Gram-Schmidt
//!# method with optional pivoting.
	return mc_orthrmxnff(m, n, a, MCLIMITS_TINYF, q, r, w, pvi);
}

MC_TARGET_FUNC int mc_rmgsmxn(const int m, const int n, const double * a, double * q, double * MC_TARGET_RESTRICT r, double * MC_TARGET_RESTRICT w, int * pvi)
{
//!# Requires a[m x n], q[m x n] and r[n x n] where 1 < n <= m.
//!# A and Q may be the same. Performing a QR decomposition using
//!# Recursive Modified Gram-Schmidt method with optional pivoting.
	return mc_orthrmxn(m, n, a, MCLIMITS_TINY, q, r, w, pvi);
}

MC_TARGET_FUNC int mc_rmgsmxnl(const int m, const int n, const long double * a, long double * q, long double * MC_TARGET_RESTRICT r, long double * MC_TARGET_RESTRICT w, int * pvi)
{
//!# Requires a[m x n], q[m x n] and r[n x n] where 1 < n <= m.
//!# A and Q may be the same. Performing a QR decomposition using
//!# Recursive Modified Gram-Schmidt method with optional pivoting.
	return mc_orthrmxnl(m, n, a, MCLIMITS_TINYL, q, r, w, pvi);
}

#endif /* !MC_RMGSMXN_H */

/* EOF */
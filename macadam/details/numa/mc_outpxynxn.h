//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_outpxynxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_outpxymxn.h>

#ifndef MC_OUTPXYNXN_H
#define MC_OUTPXYNXN_H

#pragma mark - mc_outpxynxn -

MC_TARGET_FUNC void mc_outpxynxnf(const int n, float * a, const float * x, const float * y)
{
//!# Requires a[n x n], x[n x 1] and y[n x 1].
//!# c=x*y' i.e outer product of two vectors.
	mc_outpxymxnf(n, n, a, x, y);
}

MC_TARGET_FUNC void mc_outpxynxnff(const int n, double * a, const float * x, const float * y)
{
//!# Requires a[n x n], x[n x 1] and y[n x 1].
//!# c=x*y' i.e outer product of two vectors.
	mc_outpxymxnff(n, n, a, x, y);
}

MC_TARGET_FUNC void mc_outpxynxnfd(const int n, double * a, const float * x, const double * y)
{
//!# Requires a[n x n], x[n x 1] and y[n x 1].
//!# c=x*y' i.e outer product of two vectors.
	mc_outpxymxnfd(n, n, a, x, y);
}

MC_TARGET_FUNC void mc_outpxynxndf(const int n, double * a, const double * x, const float * y)
{
//!# Requires a[n x n], x[n x 1] and y[n x 1].
//!# c=x*y' i.e outer product of two vectors.
	mc_outpxymxndf(n, n, a, x, y);
}

MC_TARGET_FUNC void mc_outpxynxn(const int n, double * a, const double * x, const double * y)
{
//!# Requires a[n x n], x[n x 1] and y[n x 1].
//!# c=x*y' i.e outer product of two vectors.
	mc_outpxymxn(n, n, a, x, y);
}

MC_TARGET_FUNC void mc_outpxynxnl(const int n, long double * a, const long double * x, const long double * y)
{
//!# Requires a[n x n], x[n x 1] and y[n x 1].
//!# c=x*y' i.e outer product of two vectors.
	mc_outpxymxnl(n, n, a, x, y);
}

#endif /* !MC_OUTPXYNXN_H */

/* EOF */
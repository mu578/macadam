//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_qrsolve3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_mulatx3x3.h>
#include <macadam/details/numa/mc_triusolve3x3.h>

#ifndef MC_QRSOLVE3X3_H
#define MC_QRSOLVE3X3_H

#pragma mark - mc_qrsolve3x3 -

MC_TARGET_FUNC int mc_qrsolve3x3f(const float q[9], const float r[9], const float p[9], const int pvi[3], const float b[3], float x[3]) 
{
//!# Solving linear system Ax=b for QR family factorization.
//!# q[m x n], r[n x n], d[n x n], p[n x n], pvi[n x 1], x[n x 1] and b[m x 1] where m=n=3.
//!# d, p and pvi can be null; for now d, p and pvi should always be, obviously null.
//!# Accepting a permutation matrix or a pivot indeces vector. Pass null accordingly.

	mc_unused(pvi);
	mc_unused(p);

	mc_mulatx3x3f(x, q, b);
	return mc_triusolve3x3f(r, x, x);
}

MC_TARGET_FUNC int mc_qrsolve3x3ff(const float q[9], const float r[9], const float p[9], const int pvi[3], const float b[3], double x[3]) 
{
//!# Solving linear system Ax=b for QR family factorization.
//!# q[m x n], r[n x n], d[n x n], p[n x n], pvi[n x 1], x[n x 1] and b[m x 1] where m=n=3.
//!# d, p and pvi can be null; for now d, p and pvi should always be, obviously null.
//!# Accepting a permutation matrix or a pivot indeces vector. Pass null accordingly.

	mc_unused(pvi);
	mc_unused(p);

	mc_mulatx3x3ff(x, q, b);
	return mc_triusolve3x3fd(r, x, x);
}

MC_TARGET_FUNC int mc_qrsolve3x3(const double q[9], const double r[9], const double p[9], const int pvi[3], const double b[3], double x[3]) 
{
//!# Solving linear system Ax=b for QR family factorization.
//!# q[m x n], r[n x n], d[n x n], p[n x n], pvi[n x 1], x[n x 1] and b[m x 1] where m=n=3.
//!# d, p and pvi can be null; for now d, p and pvi should always be, obviously null.
//!# Accepting a permutation matrix or a pivot indeces vector. Pass null accordingly.

	mc_unused(pvi);
	mc_unused(p);

	mc_mulatx3x3(x, q, b);
	return mc_triusolve3x3(r, x, x);
}

MC_TARGET_FUNC int mc_qrsolve3x3l(const long double q[9], const long double r[9], const long double p[9], const int pvi[3], const long double b[3], long double x[3]) 
{
//!# Solving linear system Ax=b for QR family factorization.
//!# q[m x n], r[n x n], d[n x n], p[n x n], pvi[n x 1], x[n x 1] and b[m x 1] where m=n=3.
//!# d, p and pvi can be null; for now d, p and pvi should always be, obviously null.
//!# Accepting a permutation matrix or a pivot indeces vector. Pass null accordingly.

	mc_unused(pvi);
	mc_unused(p);

	mc_mulatx3x3l(x, q, b);
	return mc_triusolve3x3l(r, x, x);
}

#endif /* !MC_QRSOLVE3X3_H */

/* EOF */
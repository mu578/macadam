//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_muleabmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_MULEABMXN_H
#define MC_MULEABMXN_H

#pragma mark - mc_muleabmxn -

MC_TARGET_FUNC void mc_muleabmxnf(const int m, const int n, float * MC_TARGET_RESTRICT c, const float * a, const float * b)
{
//!# Requires c[m x n], a[m x n] and b[m x n].
//!# c=a.*b i.e Hadamard product.
	int i = 0, j;
	float v, e;
	for (; i < m; i++) {
		for (j = 0; j < n; j++) {
			mc_twoproductf(a[(n * i) + j], b[(n * i) + j], &v, &e);
			c[(n * i) + j] = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_muleabmxnff(const int m, const int n, double * MC_TARGET_RESTRICT c, const float * a, const float * b)
{
//!# Requires c[m x n], a[m x n] and b[m x n].
//!# c=a.*b i.e Hadamard product.
	int i = 0, j;
	double v, e;
	for (; i < m; i++) {
		for (j = 0; j < n; j++) {
			mc_twoproduct(mc_cast(double, a[(n * i) + j]), mc_cast(double, b[(n * i) + j]), &v, &e);
			c[(n * i) + j] = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_muleabmxnfd(const int m, const int n, double * MC_TARGET_RESTRICT c, const float * a, const double * b)
{
//!# Requires c[m x n], a[m x n] and b[m x n].
//!# c=a.*b i.e Hadamard product.
	int i = 0, j;
	double v, e;
	for (; i < m; i++) {
		for (j = 0; j < n; j++) {
			mc_twoproduct(a[(n * i) + j], mc_cast(double, b[(n * i) + j]), &v, &e);
			c[(n * i) + j] = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_muleabmxndf(const int m, const int n, double * MC_TARGET_RESTRICT c, const double * a, const float * b)
{
//!# Requires c[m x n], a[m x n] and b[m x n].
//!# c=a.*b i.e Hadamard product.
	int i = 0, j;
	double v, e;
	for (; i < m; i++) {
		for (j = 0; j < n; j++) {
			mc_twoproduct(a[(n * i) + j], mc_cast(double, b[(n * i) + j]), &v, &e);
			c[(n * i) + j] = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_muleabmxn(const int m, const int n, double * MC_TARGET_RESTRICT c, const double * a, const double * b)
{
//!# Requires c[m x n], a[m x n] and b[m x n].
//!# c=a.*b i.e Hadamard product.
	int i = 0, j;
	double v, e;
	for (; i < m; i++) {
		for (j = 0; j < n; j++) {
			mc_twoproduct(a[(n * i) + j], b[(n * i) + j], &v, &e);
			c[(n * i) + j] = v + e;
		}
	}
}

MC_TARGET_FUNC void mc_muleabmxnl(const int m, const int n, long double * MC_TARGET_RESTRICT c, const long double * a, const long double * b)
{
//!# Requires c[m x n], a[m x n] and b[m x n].
//!# c=a.*b i.e Hadamard product.
	int i = 0, j;
	long double v, e;
	for (; i < m; i++) {
		for (j = 0; j < n; j++) {
			mc_twoproductl(a[(n * i) + j], b[(n * i) + j], &v, &e);
			c[(n * i) + j] = v + e;
		}
	}
}

#endif /* !MC_MULEABMXN_H */

/* EOF */
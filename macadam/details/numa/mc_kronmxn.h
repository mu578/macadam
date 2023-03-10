//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_kronmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_KRONMXN_H
#define MC_KRONMXN_H

#pragma mark - mc_kronmxn -

MC_TARGET_FUNC void mc_kronmxnf(const int m, const int n, const int p, const int q, float * MC_TARGET_RESTRICT c, const float * a, const float * b)
{
//!# Requires c[(p * m) * (q * n)], a[m x n] and b[p x q].
//!# Computing the Kronecker product of two matrices.
	int i = 0, j, k, l;	
	for(; i < m; i++){
		for(j = 0; j < n; j++){
			for(k = 0; k < p; k++) {
				for(l = 0; l < q; l++){
					c[((q * n) * ((i * p) + k)) + ((j * q) + l)] = a[(n * i) + j] * b[(q * k) + l];
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_kronmxnff(const int m, const int n, const int p, const int q, double * MC_TARGET_RESTRICT c, const float * a, const float * b)
{
//!# Requires c[(p * m) * (q * n)], a[m x n] and b[p x q].
//!# Computing the Kronecker product of two matrices.
	int i = 0, j, k, l;	
	for(; i < m; i++){
		for(j = 0; j < n; j++){
			for(k = 0; k < p; k++) {
				for(l = 0; l < q; l++){
					c[((q * n) * ((i * p) + k)) + ((j * q) + l)] = mc_cast(double, a[(n * i) + j]) * mc_cast(double, b[(q * k) + l]);
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_kronmxnfd(const int m, const int n, const int p, const int q, double * MC_TARGET_RESTRICT c, const float * a, const double * b)
{
//!# Requires c[(p * m) * (q * n)], a[m x n] and b[p x q].
//!# Computing the Kronecker product of two matrices.
	int i = 0, j, k, l;	
	for(; i < m; i++){
		for(j = 0; j < n; j++){
			for(k = 0; k < p; k++) {
				for(l = 0; l < q; l++){
					c[((q * n) * ((i * p) + k)) + ((j * q) + l)] = mc_cast(double, a[(n * i) + j]) * b[(q * k) + l];
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_kronmxndf(const int m, const int n, const int p, const int q, double * MC_TARGET_RESTRICT c, const double * a, const float * b)
{
//!# Requires c[(p * m) * (q * n)], a[m x n] and b[p x q].
//!# Computing the Kronecker product of two matrices.
	int i = 0, j, k, l;	
	for(; i < m; i++){
		for(j = 0; j < n; j++){
			for(k = 0; k < p; k++) {
				for(l = 0; l < q; l++){
					c[((q * n) * ((i * p) + k)) + ((j * q) + l)] = a[(n * i) + j] * mc_cast(double, b[(q * k) + l]);
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_kronmxn(const int m, const int n, const int p, const int q, double * MC_TARGET_RESTRICT c, const double * a, const double * b)
{
//!# Requires c[(p * m) * (q * n)], a[m x n] and b[p x q].
//!# Computing the Kronecker product of two matrices.
	int i = 0, j, k, l;	
	for(; i < m; i++){
		for(j = 0; j < n; j++){
			for(k = 0; k < p; k++) {
				for(l = 0; l < q; l++){
					c[((q * n) * ((i * p) + k)) + ((j * q) + l)] = a[(n * i) + j] * b[(q * k) + l];
				}
			}
		}
	}
}

MC_TARGET_FUNC void mc_kronmxnl(const int m, const int n, const int p, const int q, long double * MC_TARGET_RESTRICT c, const long double * a, const long double * b)
{
//!# Requires c[(p * m) * (q * n)], a[m x n] and b[p x q].
//!# Computing the Kronecker product of two matrices.
	int i = 0, j, k, l;	
	for(; i < m; i++){
		for(j = 0; j < n; j++){
			for(k = 0; k < p; k++) {
				for(l = 0; l < q; l++){
					c[((q * n) * ((i * p) + k)) + ((j * q) + l)] = a[(n * i) + j] * b[(q * k) + l];
				}
			}
		}
	}
}

#endif /* !MC_KRONMXN_H */

/* EOF */
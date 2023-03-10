//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_diag1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_zerosnxn.h>

#ifndef MC_DIAG1XN_H
#define MC_DIAG1XN_H

#pragma mark - mc_diag1xn -

MC_TARGET_FUNC void mc_diag1xnf(const int n, float * MC_TARGET_RESTRICT a, float * MC_TARGET_RESTRICT d, const int k, const int f)
{
//!# Requires a[(n + |k|) x (n + |k|)] and d[1 x n].
//!# k=0: d elements are placed on the main diagonal.
//!# k>0: d elements are placed on the +k-th superdiagonal.
//!# k<0: d elements are placed on the -k-th subdiagonal.
//!# f=0: set k-th diagonal to d elements and zeroing other elements.
//!# f=1: only set k-th diagonal to d elements.
//!# f=2: copy k-th diagonal into d.
	int i = 0, m;
	if (k > 0 ) {
		m = n + k;
		if (f != 1) { 
			mc_zerosnxnf(m, a);
		}
		if (f == 2) {
			for (; i < n; i++) {
				d[i] = a[(m * i) + (i + k)];
			}
		} else {
			for (; i < n; i++) {
				a[(m * i) + (i + k)] = d[i];
			}
		}
	} else if (k < 0 ) {
		m = n - k;
		if (f != 1) { 
			mc_zerosnxnf(m, a);
		}
		if (f == 2) {
			for (; i < n; i++) {
				d[i] = a[(m * (i - k)) + i];
			}
		} else {
			for (; i < n; i++) {
				a[(m * (i - k)) + i] = d[i];
			}
		}
	} else {
		m = n;
		if (f != 1) { 
			mc_zerosnxnf(m, a);
		}
		if (f == 2) {
			for (; i < n; i++) {
				d[i] = a[(m * i) + i];
			}
		} else {
			for (; i < n; i++) {
				a[(m * i) + i] = d[i];
			}
		}
	}
}

MC_TARGET_FUNC void mc_diag1xnff(const int n, double * MC_TARGET_RESTRICT a, float * MC_TARGET_RESTRICT d, const int k, const int f)
{
//!# Requires a[(n + |k|) x (n + |k|)] and d[1 x n].
//!# k=0: d elements are placed on the main diagonal.
//!# k>0: d elements are placed on the +k-th superdiagonal.
//!# k<0: d elements are placed on the -k-th subdiagonal.
//!# f=0: set k-th diagonal to d elements and zeroing other elements.
//!# f=1: only set k-th diagonal to d elements.
//!# f=2: copy k-th diagonal into d.
	int i = 0, m;
	if (k > 0 ) {
		m = n + k;
		if (f != 1) { 
			mc_zerosnxn(m, a);
		}
		if (f == 2) {
			for (; i < n; i++) {
				d[i] = mc_cast(float, a[(m * i) + (i + k)]);
			}
		} else {
			for (; i < n; i++) {
				a[(m * i) + (i + k)] = mc_cast(double, d[i]);
			}
		}
	} else if (k < 0 ) {
		m = n - k;
		if (f != 1) { 
			mc_zerosnxn(m, a);
		}
		if (f == 2) {
			for (; i < n; i++) {
				d[i] = mc_cast(float, a[(m * (i - k)) + i]);
			}
		} else {
			for (; i < n; i++) {
				a[(m * (i - k)) + i] = mc_cast(double, d[i]);
			}
		}
	} else {
		m = n;
		if (f != 1) { 
			mc_zerosnxn(m, a);
		}
		if (f == 2) {
			for (; i < n; i++) {
				d[i] = mc_cast(float, a[(m * i) + i]);
			}
		} else {
			for (; i < n; i++) {
				a[(m * i) + i] = mc_cast(double, d[i]);
			}
		}
	}
}

MC_TARGET_FUNC void mc_diag1xn(const int n, double * MC_TARGET_RESTRICT a, double * MC_TARGET_RESTRICT d, const int k, const int f)
{
//!# Requires a[(n + |k|) x (n + |k|)] and d[1 x n].
//!# k=0: d elements are placed on the main diagonal.
//!# k>0: d elements are placed on the +k-th superdiagonal.
//!# k<0: d elements are placed on the -k-th subdiagonal.
//!# f=0: set k-th diagonal to d elements and zeroing other elements.
//!# f=1: only set k-th diagonal to d elements.
//!# f=2: copy k-th diagonal into d.
	int i = 0, m;
	if (k > 0 ) {
		m = n + k;
		if (f != 1) { 
			mc_zerosnxn(m, a);
		}
		if (f == 2) {
			for (; i < n; i++) {
				d[i] = a[(m * i) + (i + k)];
			}
		} else {
			for (; i < n; i++) {
				a[(m * i) + (i + k)] = d[i];
			}
		}
	} else if (k < 0 ) {
		m = n - k;
		if (f != 1) { 
			mc_zerosnxn(m, a);
		}
		if (f == 2) {
			for (; i < n; i++) {
				d[i] = a[(m * (i - k)) + i];
			}
		} else {
			for (; i < n; i++) {
				a[(m * (i - k)) + i] = d[i];
			}
		}
	} else {
		m = n;
		if (f != 1) { 
			mc_zerosnxn(m, a);
		}
		if (f == 2) {
			for (; i < n; i++) {
				d[i] = a[(m * i) + i];
			}
		} else {
			for (; i < n; i++) {
				a[(m * i) + i] = d[i];
			}
		}
	}
}

MC_TARGET_FUNC void mc_diag1xnl(const int n, long double * MC_TARGET_RESTRICT a, long double * MC_TARGET_RESTRICT d, const int k, const int f)
{
//!# Requires a[(n + |k|) x (n + |k|)] and d[1 x n].
//!# k=0: d elements are placed on the main diagonal.
//!# k>0: d elements are placed on the +k-th superdiagonal.
//!# k<0: d elements are placed on the -k-th subdiagonal.
//!# f=0: set k-th diagonal to d elements and zeroing other elements.
//!# f=1: only set k-th diagonal to d elements.
//!# f=2: copy k-th diagonal into d.
	int i = 0, m;
	if (k > 0 ) {
		m = n + k;
		if (f != 1) { 
			mc_zerosnxnl(m, a);
		}
		if (f == 2) {
			for (; i < n; i++) {
				d[i] = a[(m * i) + (i + k)];
			}
		} else {
			for (; i < n; i++) {
				a[(m * i) + (i + k)] = d[i];
			}
		}
	} else if (k < 0 ) {
		m = n - k;
		if (f != 1) { 
			mc_zerosnxnl(m, a);
		}
		if (f == 2) {
			for (; i < n; i++) {
				d[i] = a[(m * (i - k)) + i];
			}
		} else {
			for (; i < n; i++) {
				a[(m * (i - k)) + i] = d[i];
			}
		}
	} else {
		m = n;
		if (f != 1) { 
			mc_zerosnxnl(m, a);
		}
		if (f == 2) {
			for (; i < n; i++) {
				d[i] = a[(m * i) + i];
			}
		} else {
			for (; i < n; i++) {
				a[(m * i) + i] = d[i];
			}
		}
	}
}

#endif /* !MC_DIAG1XN_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sytrizenxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_copynxn.h>

#ifndef MC_SYTRIZENXN_H
#define MC_SYTRIZENXN_H

#pragma mark - mc_sytrizenxn -

MC_TARGET_FUNC void mc_sytrizenxnf(const int n, float * b, const float * a, const int f)
{
//!# f=0: copy upper-triangle to lower-triangle.
//!# f=1: copy lower-triangle to upper-triangle.
	int i = 0, j;
	if (a != b) {
		mc_copynxnf(n, b, a);
	}
	for (; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i != j) {
				const int p = f == 1 ? (n * i) + j : (n * j) + i;
				const int q = f == 1 ? (n * j) + i : (n * i) + j;
				b[p]        = b[q];
			}
		}
	}
}

MC_TARGET_FUNC void mc_sytrizenxnff(const int n, double * b, const float * a, const int f)
{
//!# f=0: copy upper-triangle to lower-triangle.
//!# f=1: copy lower-triangle to upper-triangle.
	int i = 0, j;

	mc_copynxnff(n, b, a);
	for (; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i != j) {
				const int p = f == 1 ? (n * i) + j : (n * j) + i;
				const int q = f == 1 ? (n * j) + i : (n * i) + j;
				b[p]        = b[q];
			}
		}
	}
}

MC_TARGET_FUNC void mc_sytrizenxn(const int n, double * b, const double * a, const int f)
{
//!# f=0: copy upper-triangle to lower-triangle.
//!# f=1: copy lower-triangle to upper-triangle.
	int i = 0, j;
	if (a != b) {
		mc_copynxn(n, b, a);
	}
	for (; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i != j) {
				const int p = f == 1 ? (n * i) + j : (n * j) + i;
				const int q = f == 1 ? (n * j) + i : (n * i) + j;
				b[p]        = b[q];
			}
		}
	}
}

MC_TARGET_FUNC void mc_sytrizenxnl(const int n, long double * b, const long double * a, const int f)
{
//!# f=0: copy upper-triangle to lower-triangle.
//!# f=1: copy lower-triangle to upper-triangle.
	int i = 0, j;
	if (a != b) {
		mc_copynxnl(n, b, a);
	}
	for (; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i != j) {
				const int p = f == 1 ? (n * i) + j : (n * j) + i;
				const int q = f == 1 ? (n * j) + i : (n * i) + j;
				b[p]        = b[q];
			}
		}
	}
}

#endif /* !MC_SYTRIZENXN_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_flip1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcswap.h>

#ifndef MC_FLIP1XN_H
#define MC_FLIP1XN_H

#pragma mark - mc_flip1xn -

MC_TARGET_FUNC void mc_flip1xnf(const int n, float * y, const float * x)
{
//!# Requires y[n] and x[n] where 1 < n.
//!# Y and X may be the same. Reversing order of elements
//!# in vector X and storing the result into vector Y.
	int p = 0, q = n - 1;
	float w;
	if (x != y) {
		while (q >= 0) {
			y[p] = x[q]; ++p; --q;
		}
	} else {
		while (p < q) {
			mcswap_var(w, y[p], y[q]); ++p; --q;
		}
	}
}

MC_TARGET_FUNC void mc_flip1xnff(const int n, double * y, const float * x)
{
//!# Requires y[n] and x[n] where 1 < n.
//!# Y and X may be the same. Reversing order of elements
//!# in vector X and storing the result into vector Y.
	int p = 0, q = n - 1;
	while (q >= 0) {
		y[p] = mc_cast(double, x[q]); ++p; --q;
	}
}

MC_TARGET_FUNC void mc_flip1xn(const int n, double * y, const double * x)
{
//!# Requires y[n] and x[n] where 1 < n.
//!# Y and X may be the same. Reversing order of elements
//!# in vector X and storing the result into vector Y.
	int p = 0, q = n - 1;
	double w;
	if (x != y) {
		while (q >= 0) {
			y[p] = x[q]; ++p; --q;
		}
	} else {
		while (p < q) {
			mcswap_var(w, y[p], y[q]); ++p; --q;
		}
	}
}

MC_TARGET_FUNC void mc_flip1xnl(const int n, long double * y, const long double * x)
{
//!# Requires y[n] and x[n] where 1 < n.
//!# Y and X may be the same. Reversing order of elements
//!# in vector X and storing the result into vector Y.
	int p = 0, q = n - 1;
	long double w;
	if (x != y) {
		while (q >= 0) {
			y[p] = x[q]; ++p; --q;
		}
	} else {
		while (p < q) {
			mcswap_var(w, y[p], y[q]); ++p; --q;
		}
	}
}

#endif /* !MC_FLIP1XN_H */

/* EOF */
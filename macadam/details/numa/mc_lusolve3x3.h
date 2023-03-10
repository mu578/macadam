//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lusolve3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcswap.h>

#ifndef MC_LUSOLVE3X3_H
#define MC_LUSOLVE3X3_H

#pragma mark - mc_lusolve3x3 -

MC_TARGET_FUNC int mc_lusolve3x3f(const float l[9], const float u[9], const float d[9], const float p[9], const int pvi[3], const float b[3], float x[3]) 
{
//!# Solving linear system Ax=b for LU family factorization.
//!# lu[m x n], d[m x n], p[m x n], pvi[m x 1], x[m x 1] and b[m x 1] where m=n=3.
//!# d, p and pvi can be null; for now d should always be, obviously null.
//!# Accepting a permutation matrix or a pivot indeces vector. Pass null accordingly.
	const int wantpvi = mc_nonnullptr(pvi);
	const int wantp   = mc_nonnullptr(p);

	int i;
	float w;

	mc_unused(d);

//!# Computing x=b[pvi[i]] according to permutation matrix or pivot indeces vector.
	if (wantpvi) {
		x[0] = b[pvi[0]];
		x[1] = b[pvi[1]];
		x[2] = b[pvi[2]];
	} else if (wantp) {
		i    = p[0] == 1.0f ? 0 : p[1] == 1.0f ? 1 : p[2] == 1.0f ? 2 : -1;
		if (i < 0) {
			return -1;
		}
		x[0] = b[i];
		i    = p[3] == 1.0f ? 0 : p[4] == 1.0f ? 1 : p[6] == 1.0f ? 2 : -1;
		if (i < 0) {
			return -1;
		}
		x[1] = b[i];
		i    = p[6] == 1.0f ? 0 : p[7] == 1.0f ? 1 : p[8] == 1.0f ? 2 : -1;
		if (i < 0) {
			return -1;
		}
		x[2] = b[i];
	} else {
		x[0] = b[0];
		x[1] = b[1];
		x[2] = b[2];
	}

//!# Solving L*y=b.
	x[1] = x[1] - (x[0] * l[3]);
	x[2] = x[2] - (x[0] * l[6]);
	x[2] = x[2] - (x[1] * l[7]);

//!# Solving U*x=y.
	w = u[8];
	if (w != 0.0f) {
		x[2] = x[2] / w;
	}
	x[0] = x[0] - (x[2] * u[2]);
	x[1] = x[1] - (x[2] * u[5]);

	w = u[4];
	if (w != 0.0f) {
		x[1] = x[1] / w;
	}
	x[0] = x[0] - (x[1] * u[1]);

	w = u[0];
	if (w != 0.0f) {
		x[0] = x[0] / w;
	}

	return 0;
}

MC_TARGET_FUNC int mc_lusolve3x3ff(const float l[9], const float u[9], const float d[9], const float p[9], const int pvi[3], const float b[3], double x[3]) 
{
//!# Solving linear system Ax=b for LU family factorization.
//!# lu[m x n], d[m x n], p[m x n], pvi[m x 1], x[m x 1] and b[m x 1] where m=n=3.
//!# d, p and pvi can be null; for now d should always be, obviously null.
//!# Accepting a permutation matrix or a pivot indeces vector. Pass null accordingly.
	const int wantpvi = mc_nonnullptr(pvi);
	const int wantp   = mc_nonnullptr(p);

	int i;
	double w;

	mc_unused(d);

//!# Computing x=b[pvi[i]] according to permutation matrix or pivot indeces vector.
	if (wantpvi) {
		x[0] = b[pvi[0]];
		x[1] = b[pvi[1]];
		x[2] = b[pvi[2]];
	} else if (wantp) {
		i    = p[0] == 1.0f ? 0 : p[1] == 1.0f ? 1 : p[2] == 1.0f ? 2 : -1;
		if (i < 0) {
			return -1;
		}
		x[0] = mc_cast(double, b[i]);
		i    = p[3] == 1.0f ? 0 : p[4] == 1.0f ? 1 : p[6] == 1.0f ? 2 : -1;
		if (i < 0) {
			return -1;
		}
		x[1] = mc_cast(double, b[i]);
		i    = p[6] == 1.0f ? 0 : p[7] == 1.0f ? 1 : p[8] == 1.0f ? 2 : -1;
		if (i < 0) {
			return -1;
		}
		x[2] = mc_cast(double, b[i]);
	} else {
		x[0] = mc_cast(double, b[0]);
		x[1] = mc_cast(double, b[1]);
		x[2] = mc_cast(double, b[2]);
	}

//!# Solving L*y=b.
	x[1] = x[1] - (x[0] * mc_cast(double, l[3]));
	x[2] = x[2] - (x[0] * mc_cast(double, l[6]));
	x[2] = x[2] - (x[1] * mc_cast(double, l[7]));

//!# Solving U*x=y.
	w = mc_cast(double, u[8]);
	if (w != 0.0) {
		x[2] = x[2] / w;
	}
	x[0] = x[0] - (x[2] * mc_cast(double, u[2]));
	x[1] = x[1] - (x[2] * mc_cast(double, u[5]));

	w = mc_cast(double, u[4]);
	if (w != 0.0) {
		x[1] = x[1] / w;
	}
	x[0] = x[0] - (x[1] * mc_cast(double, u[1]));

	w = mc_cast(double, u[0]);
	if (w != 0.0) {
		x[0] = x[0] / w;
	}

	return 0;
}

MC_TARGET_FUNC int mc_lusolve3x3(const double l[9], const double u[9], const double d[9], const double p[9], const int pvi[3], const double b[3], double x[3]) 
{
//!# Solving linear system Ax=b for LU family factorization.
//!# lu[m x n], d[m x n], p[m x n], pvi[m x 1], x[m x 1] and b[m x 1] where m=n=3.
//!# d, p and pvi can be null; for now d should always be, obviously null.
//!# Accepting a permutation matrix or a pivot indeces vector. Pass null accordingly.
	const int wantpvi = mc_nonnullptr(pvi);
	const int wantp   = mc_nonnullptr(p);

	int i;
	double w;

	mc_unused(d);

//!# Computing x=b[pvi[i]] according to permutation matrix or pivot indeces vector.
	if (wantpvi) {
		x[0] = b[pvi[0]];
		x[1] = b[pvi[1]];
		x[2] = b[pvi[2]];
	} else if (wantp) {
		i    = p[0] == 1.0 ? 0 : p[1] == 1.0 ? 1 : p[2] == 1.0 ? 2 : -1;
		if (i < 0) {
			return -1;
		}
		x[0] = b[i];
		i    = p[3] == 1.0 ? 0 : p[4] == 1.0 ? 1 : p[6] == 1.0 ? 2 : -1;
		if (i < 0) {
			return -1;
		}
		x[1] = b[i];
		i    = p[6] == 1.0 ? 0 : p[7] == 1.0 ? 1 : p[8] == 1.0 ? 2 : -1;
		if (i < 0) {
			return -1;
		}
		x[2] = b[i];
	} else {
		x[0] = b[0];
		x[1] = b[1];
		x[2] = b[2];
	}

//!# Solving L*y=b.
	x[1] = x[1] - (x[0] * l[3]);
	x[2] = x[2] - (x[0] * l[6]);
	x[2] = x[2] - (x[1] * l[7]);

//!# Solving U*x=y.
	w = u[8];
	if (w != 0.0) {
		x[2] = x[2] / w;
	}
	x[0] = x[0] - (x[2] * u[2]);
	x[1] = x[1] - (x[2] * u[5]);

	w = u[4];
	if (w != 0.0) {
		x[1] = x[1] / w;
	}
	x[0] = x[0] - (x[1] * u[1]);

	w = u[0];
	if (w != 0.0) {
		x[0] = x[0] / w;
	}

	return 0;
}

MC_TARGET_FUNC int mc_lusolve3x3l(const long double l[9], const long double u[9], const long double d[9], const long double p[9], const int pvi[3], const long double b[3], long double x[3]) 
{
//!# Solving linear system Ax=b for LU family factorization.
//!# lu[m x n], d[m x n], p[m x n], pvi[m x 1], x[m x 1] and b[m x 1] where m=n=3.
//!# d, p and pvi can be null; for now d should always be, obviously null.
//!# Accepting a permutation matrix or a pivot indeces vector. Pass null accordingly.
	const int wantpvi = mc_nonnullptr(pvi);
	const int wantp   = mc_nonnullptr(p);

	int i;
	long double w;

	mc_unused(d);

//!# Computing x=b[pvi[i]] according to permutation matrix or pivot indeces vector.
	if (wantpvi) {
		x[0] = b[pvi[0]];
		x[1] = b[pvi[1]];
		x[2] = b[pvi[2]];
	} else if (wantp) {
		i    = p[0] == 1.0L ? 0 : p[1] == 1.0L ? 1 : p[2] == 1.0L ? 2 : -1;
		if (i < 0) {
			return -1;
		}
		x[0] = b[i];
		i    = p[3] == 1.0L ? 0 : p[4] == 1.0L ? 1 : p[6] == 1.0L ? 2 : -1;
		if (i < 0) {
			return -1;
		}
		x[1] = b[i];
		i    = p[6] == 1.0L ? 0 : p[7] == 1.0L ? 1 : p[8] == 1.0L ? 2 : -1;
		if (i < 0) {
			return -1;
		}
		x[2] = b[i];
	} else {
		x[0] = b[0];
		x[1] = b[1];
		x[2] = b[2];
	}

//!# Solving L*y=b.
	x[1] = x[1] - (x[0] * l[3]);
	x[2] = x[2] - (x[0] * l[6]);
	x[2] = x[2] - (x[1] * l[7]);

//!# Solving U*x=y.
	w = u[8];
	if (w != 0.0L) {
		x[2] = x[2] / w;
	}
	x[0] = x[0] - (x[2] * u[2]);
	x[1] = x[1] - (x[2] * u[5]);

	w = u[4];
	if (w != 0.0L) {
		x[1] = x[1] / w;
	}
	x[0] = x[0] - (x[1] * u[1]);

	w = u[0];
	if (w != 0.0L) {
		x[0] = x[0] / w;
	}

	return 0;
}

#endif /* !MC_LUSOLVE3X3_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_triusolve3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_zeros1x3.h>

#ifndef MC_TRIUSOLVE3X3_H
#define MC_TRIUSOLVE3X3_H

#pragma mark - MC_TRIUSOLVE3X3_H -

MC_TARGET_FUNC int mc_triusolve3x3f(const float u[9], const float b[3], float x[3])
{
//!# B and X may be the same. Solving the non-singular upper
//!# triangular system Ux=b, where u[3 x 3], b[3 x 1], and x[1 x 3].
	float w;
	if (x != b) {
		x[0] = b[0];
		x[1] = b[1];
		x[2] = b[2];
	}
	w = u[8];
	if (w != 0.0f) {
		x[2] = x[2] / w;
	}
	x[0] = x[0] - x[2] * u[2];
	x[1] = x[1] - x[2] * u[5];

	w = u[4];
	if (w != 0.0f) {
		x[1] = x[1] / w;
	}
	x[0] = x[0] - x[1] * u[1];

	w = u[0];
	if (w != 0.0f) {
		x[0] = x[0] / w;
	}

	return 0;
}

MC_TARGET_FUNC int mc_triusolve3x3ff(const float u[9], const float b[3], double x[3])
{
//!# B and X may be the same. Solving the non-singular upper
//!# triangular system Ux=b, where u[3 x 3], b[3 x 1], and x[1 x 3].
	double w;

	x[0] = mc_cast(double, b[0]);
	x[1] = mc_cast(double, b[1]);
	x[2] = mc_cast(double, b[2]);

	w = mc_cast(double, u[8]);
	if (w != 0.0) {
		x[2] = x[2] / w;
	}
	x[0] = x[0] - x[2] * mc_cast(double, u[2]);
	x[1] = x[1] - x[2] * mc_cast(double, u[5]);

	w = mc_cast(double, u[4]);
	if (w != 0.0) {
		x[1] = x[1] / w;
	}
	x[0] = x[0] - x[1] * mc_cast(double, u[1]);

	w = mc_cast(double, u[0]);
	if (w != 0.0) {
		x[0] = x[0] / w;
	}

	return 0;
}

MC_TARGET_FUNC int mc_triusolve3x3fd(const float u[9], const double b[3], double x[3])
{
//!# B and X may be the same. Solving the non-singular upper
//!# triangular system Ux=b, where u[3 x 3], b[3 x 1], and x[1 x 3].
	double w;
	if (x != b) {
		x[0] = b[0];
		x[1] = b[1];
		x[2] = b[2];
	}
	w = mc_cast(double, u[8]);
	if (w != 0.0) {
		x[2] = x[2] / w;
	}
	x[0] = x[0] - x[2] * mc_cast(double, u[2]);
	x[1] = x[1] - x[2] * mc_cast(double, u[5]);

	w = mc_cast(double, u[4]);
	if (w != 0.0) {
		x[1] = x[1] / w;
	}
	x[0] = x[0] - x[1] * mc_cast(double, u[1]);

	w = mc_cast(double, u[0]);
	if (w != 0.0) {
		x[0] = x[0] / w;
	}

	return 0;
}

MC_TARGET_FUNC int mc_triusolve3x3(const double u[9], const double b[3], double x[3])
{
//!# B and X may be the same. Solving the non-singular upper
//!# triangular system Ux=b, where u[3 x 3], b[3 x 1], and x[1 x 3].
	double w;
	if (x != b) {
		x[0] = b[0];
		x[1] = b[1];
		x[2] = b[2];
	}
	w = u[8];
	if (w != 0.0) {
		x[2] = x[2] / w;
	}
	x[0] = x[0] - x[2] * u[2];
	x[1] = x[1] - x[2] * u[5];

	w = u[4];
	if (w != 0.0) {
		x[1] = x[1] / w;
	}
	x[0] = x[0] - x[1] * u[1];

	w = u[0];
	if (w != 0.0) {
		x[0] = x[0] / w;
	}

	return 0;
}

MC_TARGET_FUNC int mc_triusolve3x3l(const long double u[9], const long double b[3], long double x[3])
{
//!# B and X may be the same. Solving the non-singular upper
//!# triangular system Ux=b, where u[3 x 3], b[3 x 1], and x[1 x 3].
	long double w;
	if (x != b) {
		x[0] = b[0];
		x[1] = b[1];
		x[2] = b[2];
	}
	w = u[8];
	if (w != 0.0L) {
		x[2] = x[2] / w;
	}
	x[0] = x[0] - x[2] * u[2];
	x[1] = x[1] - x[2] * u[5];

	w = u[4];
	if (w != 0.0L) {
		x[1] = x[1] / w;
	}
	x[0] = x[0] - x[1] * u[1];

	w = u[0];
	if (w != 0.0L) {
		x[0] = x[0] / w;
	}

	return 0;
}

#endif /* !MC_TRIUSOLVE3X3_H */

/* EOF */
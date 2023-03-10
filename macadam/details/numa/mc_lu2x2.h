//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lu2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_eye2x2.h>
#include <macadam/details/numa/mc_zeros2x2.h>

#ifndef MC_LU2X2_H
#define MC_LU2X2_H

#pragma mark - mc_lu2x2 -

MC_TARGET_FUNC int mc_lu2x2f(const float a[4], float l[4], float u[4])
{
//!# A and L may be the same. Returns A=L(DU) as per Doolittle's method.
	float a11 = a[0], a12 = a[1];
	float a21 = a[2], a22 = a[3];

	mc_eye2x2f(l, 0);
	mc_zeros2x2f(u);

	u[0] = a11;
	u[1] = a12;
	if (u[0] != 0.0f) {
		l[2] = a21 / u[0];
	}
	u[3] = a22 - (l[2] * u[1]);

	return 0;
}

MC_TARGET_FUNC int mc_lu2x2ff(const float a[4], double l[4], double u[4])
{
//!# A and L may be the same. Returns A=L(DU) as per Doolittle's method.
	double a11 = mc_cast(double, a[0]), a12 = mc_cast(double, a[1]);
	double a21 = mc_cast(double, a[2]), a22 = mc_cast(double, a[3]);

	mc_eye2x2(l, 0);
	mc_zeros2x2(u);

	u[0] = a11;
	u[1] = a12;
	if (u[0] != 0.0) {
		l[2] = a21 / u[0];
	}
	u[3] = a22 - (l[2] * u[1]);

	return 0;
}

MC_TARGET_FUNC int mc_lu2x2(const double a[4], double l[4], double u[4])
{
//!# A and L may be the same. Returns A=L(DU) as per Doolittle's method.
	double a11 = a[0], a12 = a[1];
	double a21 = a[2], a22 = a[3];

	mc_eye2x2(l, 0);
	mc_zeros2x2(u);

	u[0] = a11;
	u[1] = a12;
	if (u[0] != 0.0) {
		l[2] = a21 / u[0];
	}
	u[3] = a22 - (l[2] * u[1]);

	return 0;
}

MC_TARGET_FUNC int mc_lu2x2l(const long double a[4], long double l[4], long double u[4])
{
//!# A and L may be the same. Returns A=L(DU) as per Doolittle's method.
	long double a11 = a[0], a12 = a[1];
	long double a21 = a[2], a22 = a[3];

	mc_eye2x2l(l, 0);
	mc_zeros2x2l(u);

	u[0] = a11;
	u[1] = a12;
	if (u[0] != 0.0L) {
		l[2] = a21 / u[0];
	}
	u[3] = a22 - (l[2] * u[1]);

	return 0;
}

#endif /* !MC_LU2X2_H */

/* EOF */
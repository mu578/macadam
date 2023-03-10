//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ldu3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_eye3x3.h>
#include <macadam/details/numa/mc_lu3x3.h>

#ifndef MC_LDU3X3_H
#define MC_LDU3X3_H

#pragma mark - mc_ldu3x3 -

MC_TARGET_FUNC int mc_ldu3x3f(const float a[9], float l[9], float d[9], float u[9])
{
//!# A and L may be the same. Using a closed-form expression.
	mc_eye3x3f(d, 0);
	if (0 == mc_lu3x3f(a, l, u)) {
		if (u[0] != 0.0f && u[4] != 0.0f) {
//!# Updating diagonal.
			d[0] = u[0]; d[4] = u[4]; d[8] = u[8];
//!# Updating U with REF of upper triangle.
			u[0] = 1.0f; u[1] = u[1] / d[0]; u[2] = u[2] / d[0];
			u[3] = 0.0f; u[4] = 1.0f;        u[5] = u[5] / d[4];
			u[6] = 0.0f; u[7] = 0.0f;        u[8] = 1.0f;
			return 0;
		}
		mc_eye3x3f(l, 0);
		mc_eye3x3f(u, 0);
	}
	return -1;
}

MC_TARGET_FUNC int mc_ldu3x3ff(const float a[9], double l[9], double d[9], double u[9])
{
//!# Using a closed-form expression.
	mc_eye3x3(d, 0);
	if (0 == mc_lu3x3ff(a, l, u)) {
		if (u[0] != 0.0 && u[4] != 0.0) {
//!# Updating diagonal.
			d[0] = u[0]; d[4] = u[4]; d[8] = u[8];
//!# Updating U with REF of upper triangle.
			u[0] = 1.0; u[1] = u[1] / d[0]; u[2] = u[2] / d[0];
			u[3] = 0.0; u[4] = 1.0;         u[5] = u[5] / d[4];
			u[6] = 0.0; u[7] = 0.0;         u[8] = 1.0;
			return 0;
		}
		mc_eye3x3(l, 0);
		mc_eye3x3(u, 0);
	}
	return -1;
}

MC_TARGET_FUNC int mc_ldu3x3(const double a[9], double l[9], double d[9], double u[9])
{
//!# A and L may be the same. Using a closed-form expression.
	mc_eye3x3(d, 0);
	if (0 == mc_lu3x3(a, l, u)) {
		if (u[0] != 0.0 && u[4] != 0.0) {
//!# Updating diagonal.
			d[0] = u[0]; d[4] = u[4]; d[8] = u[8];
//!# Updating U with REF of upper triangle.
			u[0] = 1.0; u[1] = u[1] / d[0]; u[2] = u[2] / d[0];
			u[3] = 0.0; u[4] = 1.0;         u[5] = u[5] / d[4];
			u[6] = 0.0; u[7] = 0.0;         u[8] = 1.0;
			return 0;
		}
		mc_eye3x3(l, 0);
		mc_eye3x3(u, 0);
	}
	return -1;
}

MC_TARGET_FUNC int mc_ldu3x3l(const long double a[9], long double l[9], long double d[9], long double u[9])
{
//!# A and L may be the same. Using a closed-form expression.
	mc_eye3x3l(d, 0);
	if (0 == mc_lu3x3l(a, l, u)) {
		if (u[0] != 0.0L && u[4] != 0.0L) {
//!# Updating diagonal.
			d[0] = u[0]; d[4] = u[4]; d[8] = u[8];
//!# Updating U with REF of upper triangle.
			u[0] = 1.0L; u[1] = u[1] / d[0]; u[2] = u[2] / d[0];
			u[3] = 0.0L; u[4] = 1.0L;        u[5] = u[5] / d[4];
			u[6] = 0.0L; u[7] = 0.0L;        u[8] = 1.0L;
			return 0;
		}
		mc_eye3x3l(l, 0);
		mc_eye3x3l(u, 0);
	}
	return -1;
}

#endif /* !MC_LDU3X3_H */

/* EOF */
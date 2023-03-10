//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_qrgs2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_hypot2.h>
#include <macadam/details/numa/mc_eye2x2.h>
#include <macadam/details/numa/mc_zeros2x2.h>

#ifndef MC_QRGS2X2_H
#define MC_QRGS2X2_H

#pragma mark - mc_qrgs2x2 -

MC_TARGET_FUNC int mc_qrgs2x2f(const float a[4], float q[4], float r[4])
{
//!# A and Q may be the same. Using Gram-Schmidt method.
	float w;

	if (a != q) {
		q[0] = a[0]; q[1] = a[1];
		q[2] = a[2]; q[3] = a[3];
	}
	mc_zeros2x2f(r);

	w    = mc_hypot2f(q[0], q[2]);
	r[0] = w;
	if (w != 0.0f) {
		w    = 1.0f / w;
		q[0] = q[0] * w;
		q[2] = q[2] * w;
	} else {
		q[0] = 1.0f;
	}

	w    = q[0] * q[1] + q[2] * q[3];
	r[1] = w;

	q[1] = q[1] - w * q[0];
	q[3] = q[3] - w * q[2];

	w    = mc_hypot2f(q[1], q[3]);
	r[3] = w;
	if (w != 0.0f) {
		w    = 1.0f / w;
		q[1] = q[1] * w;
		q[3] = q[3] * w;
	} else {
		q[1] = 1.0f;
	}

	return 0;
}

MC_TARGET_FUNC int mc_qrgs2x2ff(const float a[4], double q[4], double r[4])
{
//!# Using Gram-Schmidt method.
	double w;

	q[0] = mc_cast(double, a[0]); q[1] =  mc_cast(double, a[1]);
	q[2] = mc_cast(double, a[2]); q[3] =  mc_cast(double, a[3]);

	mc_zeros2x2(r);

	w    = mc_hypot2(q[0], q[2]);
	r[0] = w;
	if (w != 0.0) {
		w    = 1.0 / w;
		q[0] = q[0] * w;
		q[2] = q[2] * w;
	} else {
		q[0] = 1.0;
	}

	w    = q[0] * q[1] + q[2] * q[3];
	r[1] = w;

	q[1] = q[1] - w * q[0];
	q[3] = q[3] - w * q[2];

	w    = mc_hypot2(q[1], q[3]);
	r[3] = w;
	if (w != 0.0) {
		w    = 1.0 / w;
		q[1] = q[1] * w;
		q[3] = q[3] * w;
	} else {
		q[1] = 1.0;
	}

	return 0;
}

MC_TARGET_FUNC int mc_qrgs2x2(const double a[4], double q[4], double r[4])
{
//!# A and Q may be the same. Using Gram-Schmidt method.
	double w;

	if (a != q) {
		q[0] = a[0]; q[1] = a[1];
		q[2] = a[2]; q[3] = a[3];
	}

	mc_zeros2x2(r);

	w    = mc_hypot2(q[0], q[2]);
	r[0] = w;
	if (w != 0.0) {
		w    = 1.0 / w;
		q[0] = q[0] * w;
		q[2] = q[2] * w;
	} else {
		q[0] = 1.0;
	}

	w    = q[0] * q[1] + q[2] * q[3];
	r[1] = w;

	q[1] = q[1] - w * q[0];
	q[3] = q[3] - w * q[2];

	w    = mc_hypot2(q[1], q[3]);
	r[3] = w;
	if (w != 0.0) {
		w    = 1.0 / w;
		q[1] = q[1] * w;
		q[3] = q[3] * w;
	} else {
		q[1] = 1.0;
	}

	return 0;
}

MC_TARGET_FUNC int mc_qrgs2x2l(const long double a[4], long double q[4], long double r[4])
{
//!# A and Q may be the same. Using Gram-Schmidt method.
	long double w;

	if (a != q) {
		q[0] = a[0]; q[1] = a[1];
		q[2] = a[2]; q[3] = a[3];
	}

	mc_zeros2x2l(r);

	w    = mc_hypot2l(q[0], q[2]);
	r[0] = w;
	if (w != 0.0L) {
		w    = 1.0L / w;
		q[0] = q[0] * w;
		q[2] = q[2] * w;
	} else {
		q[0] = 1.0L;
	}

	w    = q[0] * q[1] + q[2] * q[3];
	r[1] = w;

	q[1] = q[1] - w * q[0];
	q[3] = q[3] - w * q[2];

	w    = mc_hypot2l(q[1], q[3]);
	r[3] = w;
	if (w != 0.0L) {
		w    = 1.0L / w;
		q[1] = q[1] * w;
		q[3] = q[3] * w;
	} else {
		q[1] = 1.0L;
	}

	return 0;
}

#endif /* !MC_QRGS2X2_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_qrgs3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_hypot3.h>
#include <macadam/details/numa/mc_zeros3x3.h>

#ifndef MC_QRGS3X3_H
#define MC_QRGS3X3_H

#pragma mark - mc_qrgs3x3 -

MC_TARGET_FUNC int mc_qrgs3x3f(const float a[9], float q[9], float r[9])
{
//!# A and Q may be the same. QR using Modified Gram-Schmidt
//!# method, A must be fullrank.
	float w;

	if (a != q) {
		q[0] = a[0]; q[1] = a[1]; q[2] = a[2];
		q[3] = a[3]; q[4] = a[4]; q[5] = a[5];
		q[6] = a[6]; q[7] = a[7]; q[8] = a[8];
	}

	mc_zeros3x3f(r);

	w    = mc_hypot3f(q[0], q[3], q[6]);
	r[0] = w;
	if (w != 0.0f) {
		w    = 1.0f / w;
		q[0] = q[0] * w;
		q[3] = q[3] * w;
		q[6] = q[6] * w;
	} else {
		q[0] = 1.0f;
	}

	w    = q[0] * q[1] + q[3] * q[4] + q[6] * q[7];
	r[1] = w;

	q[1] = q[1] - w * q[0];
	q[4] = q[4] - w * q[3];
	q[7] = q[7] - w * q[6];

	w    = mc_hypot3f(q[1], q[4], q[7]);
	r[4] = w;
	if (w != 0.0f) {
		w    = 1.0f / w;
		q[1] = q[1] * w;
		q[4] = q[4] * w;
		q[7] = q[7] * w;
	} else {
		q[1] = 1.0f;
	}

	w    = q[0] * q[2] + q[3] * q[5] + q[6] * q[8];
	r[2] = w;

	q[2] = q[2] - w * q[0];
	q[5] = q[5] - w * q[3];
	q[8] = q[8] - w * q[6];

	w    = q[1] * q[2] + q[4] * q[5] + q[7] * q[8];
	r[5] = w;

	q[2] = q[2] - w * q[1];
	q[5] = q[5] - w * q[4];
	q[8] = q[8] - w * q[7];

	w    = mc_hypot3f(q[2], q[5], q[8]);
	r[8] = w;
	if (w != 0.0f) {
		w    = 1.0f / w;
		q[2] = q[2] * w;
		q[5] = q[5] * w;
		q[8] = q[8] * w;
	} else {
		q[2] = 1.0f;
	}

	return 0;
}

MC_TARGET_FUNC int mc_qrgs3x3ff(const float a[9], double q[9], double r[9])
{
//!# QR using Modified Gram-Schmidt method, A must be fullrank.
	double w;

	q[0] = mc_cast(double, a[0]); q[1] = mc_cast(double, a[1]); q[2] = mc_cast(double, a[2]);
	q[3] = mc_cast(double, a[3]); q[4] = mc_cast(double, a[4]); q[5] = mc_cast(double, a[5]);
	q[6] = mc_cast(double, a[6]); q[7] = mc_cast(double, a[7]); q[8] = mc_cast(double, a[8]);

	mc_zeros3x3(r);

	w    = mc_hypot3(q[0], q[3], q[6]);
	r[0] = w;
	if (w != 0.0) {
		w    = 1.0 / w;
		q[0] = q[0] * w;
		q[3] = q[3] * w;
		q[6] = q[6] * w;
	} else {
		q[0] = 1.0;
	}

	w    = q[0] * q[1] + q[3] * q[4] + q[6] * q[7];
	r[1] = w;

	q[1] = q[1] - w * q[0];
	q[4] = q[4] - w * q[3];
	q[7] = q[7] - w * q[6];

	w    = mc_hypot3(q[1], q[4], q[7]);
	r[4] = w;
	if (w != 0.0) {
		w    = 1.0 / w;
		q[1] = q[1] * w;
		q[4] = q[4] * w;
		q[7] = q[7] * w;
	} else {
		q[1] = 1.0;
	}

	w    = q[0] * q[2] + q[3] * q[5] + q[6] * q[8];
	r[2] = w;

	q[2] = q[2] - w * q[0];
	q[5] = q[5] - w * q[3];
	q[8] = q[8] - w * q[6];

	w    = q[1] * q[2] + q[4] * q[5] + q[7] * q[8];
	r[5] = w;

	q[2] = q[2] - w * q[1];
	q[5] = q[5] - w * q[4];
	q[8] = q[8] - w * q[7];

	w    = mc_hypot3(q[2], q[5], q[8]);
	r[8] = w;
	if (w != 0.0) {
		w    = 1.0 / w;
		q[2] = q[2] * w;
		q[5] = q[5] * w;
		q[8] = q[8] * w;
	} else {
		q[2] = 1.0;
	}

	return 0;
}

MC_TARGET_FUNC int mc_qrgs3x3(const double a[9], double q[9], double r[9])
{
//!# A and Q may be the same. QR using Modified Gram-Schmidt
//!# method, A must be fullrank.
	double w;

	if (a != q) {
		q[0] = a[0]; q[1] = a[1]; q[2] = a[2];
		q[3] = a[3]; q[4] = a[4]; q[5] = a[5];
		q[6] = a[6]; q[7] = a[7]; q[8] = a[8];
	}

	mc_zeros3x3(r);

	w    = mc_hypot3(q[0], q[3], q[6]);
	r[0] = w;
	if (w != 0.0) {
		w    = 1.0 / w;
		q[0] = q[0] * w;
		q[3] = q[3] * w;
		q[6] = q[6] * w;
	} else {
		q[0] = 1.0;
	}

	w    = q[0] * q[1] + q[3] * q[4] + q[6] * q[7];
	r[1] = w;

	q[1] = q[1] - w * q[0];
	q[4] = q[4] - w * q[3];
	q[7] = q[7] - w * q[6];

	w    = mc_hypot3(q[1], q[4], q[7]);
	r[4] = w;
	if (w != 0.0) {
		w    = 1.0 / w;
		q[1] = q[1] * w;
		q[4] = q[4] * w;
		q[7] = q[7] * w;
	} else {
		q[1] = 1.0;
	}

	w    = q[0] * q[2] + q[3] * q[5] + q[6] * q[8];
	r[2] = w;

	q[2] = q[2] - w * q[0];
	q[5] = q[5] - w * q[3];
	q[8] = q[8] - w * q[6];

	w    = q[1] * q[2] + q[4] * q[5] + q[7] * q[8];
	r[5] = w;

	q[2] = q[2] - w * q[1];
	q[5] = q[5] - w * q[4];
	q[8] = q[8] - w * q[7];

	w    = mc_hypot3(q[2], q[5], q[8]);
	r[8] = w;
	if (w != 0.0) {
		w    = 1.0 / w;
		q[2] = q[2] * w;
		q[5] = q[5] * w;
		q[8] = q[8] * w;
	} else {
		q[2] = 1.0;
	}

	return 0;
}

MC_TARGET_FUNC int mc_qrgs3x3l(const long double a[9], long double q[9], long double r[9])
{
//!# A and Q may be the same. QR using Modified Gram-Schmidt
//!# method, A must be fullrank.
	long double w;

	if (a != q) {
		q[0] = a[0]; q[1] = a[1]; q[2] = a[2];
		q[3] = a[3]; q[4] = a[4]; q[5] = a[5];
		q[6] = a[6]; q[7] = a[7]; q[8] = a[8];
	}

	mc_zeros3x3l(r);

	w    = mc_hypot3l(q[0], q[3], q[6]);
	r[0] = w;
	if (w != 0.0L) {
		w    = 1.0L / w;
		q[0] = q[0] * w;
		q[3] = q[3] * w;
		q[6] = q[6] * w;
	} else {
		q[0] = 1.0L;
	}

	w    = q[0] * q[1] + q[3] * q[4] + q[6] * q[7];
	r[1] = w;

	q[1] = q[1] - w * q[0];
	q[4] = q[4] - w * q[3];
	q[7] = q[7] - w * q[6];

	w    = mc_hypot3l(q[1], q[4], q[7]);
	r[4] = w;
	if (w != 0.0L) {
		w    = 1.0L / w;
		q[1] = q[1] * w;
		q[4] = q[4] * w;
		q[7] = q[7] * w;
	} else {
		q[1] = 1.0L;
	}

	w    = q[0] * q[2] + q[3] * q[5] + q[6] * q[8];
	r[2] = w;

	q[2] = q[2] - w * q[0];
	q[5] = q[5] - w * q[3];
	q[8] = q[8] - w * q[6];

	w    = q[1] * q[2] + q[4] * q[5] + q[7] * q[8];
	r[5] = w;

	q[2] = q[2] - w * q[1];
	q[5] = q[5] - w * q[4];
	q[8] = q[8] - w * q[7];

	w    = mc_hypot3l(q[2], q[5], q[8]);
	r[8] = w;
	if (w != 0.0L) {
		w    = 1.0L / w;
		q[2] = q[2] * w;
		q[5] = q[5] * w;
		q[8] = q[8] * w;
	} else {
		q[2] = 1.0L;
	}

	return 0;
}

#endif /* !MC_QRGS3X3_H */

/* EOF */
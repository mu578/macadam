//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zreig2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fmax.h>
#include <macadam/details/math/mc_fmin.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/math/mc_zadd.h>
#include <macadam/details/math/mc_zmul.h>
#include <macadam/details/math/mc_zsqrt.h>
#include <macadam/details/math/mc_zsub.h>
#include <macadam/details/numa/mc_det2x2.h>
#include <macadam/details/numa/mc_trace2x2.h>
#include <macadam/details/numa/mc_zeye2x2.h>
#include <macadam/details/numa/mc_zunit1x2.h>

#ifndef MC_ZREIG2X2_H
#define MC_ZREIG2X2_H

#pragma mark - mc_zreig2x2 -

MC_TARGET_PROC int mc_zreig2x2f(
	  const float a[4]
	, float *     e0_r
	, float *     e0_i
	, float *     e1_r
	, float *     e1_i
	, int         wantv
	, float * v0_r, float * v0_i, float * v1_r, float * v1_i
	, float * v2_r, float * v2_i, float * v3_r, float * v3_i
) {
	float t0, t1, t2;

	if (a[2] == 0.0f && a[1] == 0.0f) {
		*e0_r = mc_fminf(a[0], a[3]); *e0_i = 0.0f;
		*e1_r = mc_fmaxf(a[0], a[3]); *e1_i = 0.0f;
		if (wantv) {
			mc_zeye2x2f(
				  v0_r, v0_i, v1_r, v1_i
				, v2_r, v2_i, v3_r, v3_i
				, 0
			);
		}
		return 2;
	} else {
		const float t = mc_trace2x2f(a);
		const float d = mc_det2x2f(a);

		mc_zsqrtf(&t0, &t1, mc_raise2f(t) / 4.0f - d, 0.0f);
		t2 = t0;

		mc_zaddf(e0_r, e0_i, t * 0.5f, 0.0f, 0.0f, t1);
		mc_zsubf(e1_r, e1_i, t * 0.5f, 0.0f, 0.0f, t1);

		if (wantv) {
			if (a[2] != 0.0f && mc_fabsf(a[2]) > mc_fabsf(a[1])) {
				mc_zsubf(&t0, &t1, *e0_r, *e0_i, a[3], 0.0f);
				*v0_r = t0; *v0_i = t1;

				mc_zsubf(&t0, &t1, *e1_r, *e1_i, a[3], 0.0f);
				*v1_r = t0; *v1_i = t1;

				*v2_r = a[2]; *v2_i = 0.0f; *v3_r = a[2]; *v3_i = 0.0f;

				mc_zunit1x2f(v0_r, v0_i, v2_r, v2_i);
				mc_zunit1x2f(v1_r, v1_i, v3_r, v3_i);

			} else if (a[1] != 0.0f) {
				*v0_r = a[1]; *v0_i = 0.0f; *v1_r = a[1]; *v1_i = 0.0f;

				mc_zsubf(&t0, &t1, *e0_r, *e0_i, a[0], 0.0f);
				*v2_r = t0; *v2_i = t1;

				mc_zsubf(&t0, &t1, *e1_r, *e1_i, a[0], 0.0f);
				*v3_r = t0; *v3_i = t1;

				mc_zunit1x2f(v0_r, v0_i, v2_r, v2_i);
				mc_zunit1x2f(v1_r, v1_i, v3_r, v3_i);
			} else {
				mc_zeye2x2f(
					  v0_r, v0_i, v1_r, v1_i
					, v2_r, v2_i, v3_r, v3_i
					, 0
				);
			}
		}
		return a[2] == a[1] ? 1 : (t2 > 0.0f ? 0 : 1);
	}
	return -1;
}

MC_TARGET_PROC int mc_zreig2x2ff(
	  const float a[4]
	, double *    e0_r
	, double *    e0_i
	, double *    e1_r
	, double *    e1_i
	, int         wantv
	, double * v0_r, double * v0_i, double * v1_r, double * v1_i
	, double * v2_r, double * v2_i, double * v3_r, double * v3_i
) {
	double t0, t1, t2;

	if (a[2] == 0.0f && a[1] == 0.0f) {
		*e0_r = mc_cast(double, mc_fminf(a[0], a[3])); *e0_i = 0.0;
		*e1_r = mc_cast(double, mc_fmaxf(a[0], a[3])); *e1_i = 0.0;
		if (wantv) {
			mc_zeye2x2(
				  v0_r, v0_i, v1_r, v1_i
				, v2_r, v2_i, v3_r, v3_i
				, 0
			);
		}
		return 2;
	} else {
		const double t = mc_trace2x2ff(a);
		const double d = mc_det2x2ff(a);

		mc_zsqrt(&t0, &t1, mc_raise2(t) / 4.0 - d, 0.0);
		t2 = t0;

		mc_zadd(e0_r, e0_i, t * 0.5, 0.0, 0.0, t1);
		mc_zsub(e1_r, e1_i, t * 0.5, 0.0, 0.0, t1);

		if (wantv) {
			if (a[2] != 0.0f && mc_fabsf(a[2]) > mc_fabsf(a[1])) {
				mc_zsub(&t0, &t1, *e0_r, *e0_i, mc_cast(double, a[3]), 0.0);
				*v0_r = t0; *v0_i = t1;

				mc_zsub(&t0, &t1, *e1_r, *e1_i, mc_cast(double, a[3]), 0.0);
				*v1_r = t0; *v1_i = t1;

				*v2_r = mc_cast(double, a[2]); *v2_i = 0.0; *v3_r = mc_cast(double, a[2]); *v3_i = 0.0;

				mc_zunit1x2(v0_r, v0_i, v2_r, v2_i);
				mc_zunit1x2(v1_r, v1_i, v3_r, v3_i);

			} else if (a[1] != 0.0f) {
				*v0_r = mc_cast(double, a[1]); *v0_i = 0.0; *v1_r = mc_cast(double, a[1]); *v1_i = 0.0;

				mc_zsub(&t0, &t1, *e0_r, *e0_i, mc_cast(double, a[0]), 0.0);
				*v2_r = t0; *v2_i = t1;

				mc_zsub(&t0, &t1, *e1_r, *e1_i, mc_cast(double, a[0]), 0.0);
				*v3_r = t0; *v3_i = t1;

				mc_zunit1x2(v0_r, v0_i, v2_r, v2_i);
				mc_zunit1x2(v1_r, v1_i, v3_r, v3_i);
			} else {
				mc_zeye2x2(
					  v0_r, v0_i, v1_r, v1_i
					, v2_r, v2_i, v3_r, v3_i
					, 0
				);
			}
		}
		return a[2] == a[1] ? 1 : (t2 > 0.0 ? 0 : 1);
	}
	return -1;
}

MC_TARGET_PROC int mc_zreig2x2(
	  const double a[4]
	, double *     e0_r
	, double *     e0_i
	, double *     e1_r
	, double *     e1_i
	, int          wantv
	, double * v0_r, double * v0_i, double * v1_r, double * v1_i
	, double * v2_r, double * v2_i, double * v3_r, double * v3_i
) {
	double t0, t1, t2;

	if (a[2] == 0.0 && a[1] == 0.0) {
		*e0_r = mc_fmin(a[0], a[3]); *e0_i = 0.0;
		*e1_r = mc_fmax(a[0], a[3]); *e1_i = 0.0;
		if (wantv) {
			mc_zeye2x2(
				  v0_r, v0_i, v1_r, v1_i
				, v2_r, v2_i, v3_r, v3_i
				, 0
			);
		}
		return 2;
	} else {
		const double t = mc_trace2x2(a);
		const double d = mc_det2x2(a);

		mc_zsqrt(&t0, &t1, mc_raise2(t) / 4.0 - d, 0.0);
		t2 = t0;

		mc_zadd(e0_r, e0_i, t * 0.5, 0.0, 0.0, t1);
		mc_zsub(e1_r, e1_i, t * 0.5, 0.0, 0.0, t1);

		if (wantv) {
			if (a[2] != 0.0 && mc_fabs(a[2]) > mc_fabs(a[1])) {
				mc_zsub(&t0, &t1, *e0_r, *e0_i, a[3], 0.0);
				*v0_r = t0; *v0_i = t1;

				mc_zsub(&t0, &t1, *e1_r, *e1_i, a[3], 0.0);
				*v1_r = t0; *v1_i = t1;

				*v2_r = a[2]; *v2_i = 0.0; *v3_r = a[2]; *v3_i = 0.0;

				mc_zunit1x2(v0_r, v0_i, v2_r, v2_i);
				mc_zunit1x2(v1_r, v1_i, v3_r, v3_i);

			} else if (a[1] != 0.0) {
				*v0_r = a[1]; *v0_i = 0.0; *v1_r = a[1]; *v1_i = 0.0;

				mc_zsub(&t0, &t1, *e0_r, *e0_i, a[0], 0.0);
				*v2_r = t0; *v2_i = t1;

				mc_zsub(&t0, &t1, *e1_r, *e1_i, a[0], 0.0);
				*v3_r = t0; *v3_i = t1;

				mc_zunit1x2(v0_r, v0_i, v2_r, v2_i);
				mc_zunit1x2(v1_r, v1_i, v3_r, v3_i);
			} else {
				mc_zeye2x2(
					  v0_r, v0_i, v1_r, v1_i
					, v2_r, v2_i, v3_r, v3_i
					, 0
				);
			}
		}
		return a[2] == a[1] ? 1 : (t2 > 0.0 ? 0 : 1);
	}
	return -1;
}

MC_TARGET_PROC int mc_zreig2x2l(
	  const long double a[4]
	, long double *     e0_r
	, long double *     e0_i
	, long double *     e1_r
	, long double *     e1_i
	, int               wantv
	, long double * v0_r, long double * v0_i, long double * v1_r, long double * v1_i
	, long double * v2_r, long double * v2_i, long double * v3_r, long double * v3_i
) {
	long double t0, t1, t2;

	if (a[2] == 0.0L && a[1] == 0.0L) {
		*e0_r = mc_fminl(a[0], a[3]); *e0_i = 0.0L;
		*e1_r = mc_fmaxl(a[0], a[3]); *e1_i = 0.0L;
		if (wantv) {
			mc_zeye2x2l(
				  v0_r, v0_i, v1_r, v1_i
				, v2_r, v2_i, v3_r, v3_i
				, 0
			);
		}
		return 2;
	} else {
		const long double t = mc_trace2x2l(a);
		const long double d = mc_det2x2l(a);

		mc_zsqrtl(&t0, &t1, mc_raise2l(t) / 4.0L - d, 0.0L);
		t2 = t0;

		mc_zaddl(e0_r, e0_i, t * 0.5L, 0.0L, 0.0L, t1);
		mc_zsubl(e1_r, e1_i, t * 0.5L, 0.0L, 0.0L, t1);

		if (wantv) {
			if (a[2] != 0.0L && mc_fabsl(a[2]) > mc_fabsl(a[1])) {
				mc_zsubl(&t0, &t1, *e0_r, *e0_i, a[3], 0.0L);
				*v0_r = t0; *v0_i = t1;

				mc_zsubl(&t0, &t1, *e1_r, *e1_i, a[3], 0.0L);
				*v1_r = t0; *v1_i = t1;

				*v2_r = a[2]; *v2_i = 0.0L; *v3_r = a[2]; *v3_i = 0.0L;

				mc_zunit1x2l(v0_r, v0_i, v2_r, v2_i);
				mc_zunit1x2l(v1_r, v1_i, v3_r, v3_i);

			} else if (a[1] != 0.0L) {
				*v0_r = a[1]; *v0_i = 0.0L; *v1_r = a[1]; *v1_i = 0.0L;

				mc_zsubl(&t0, &t1, *e0_r, *e0_i, a[0], 0.0L);
				*v2_r = t0; *v2_i = t1;

				mc_zsubl(&t0, &t1, *e1_r, *e1_i, a[0], 0.0L);
				*v3_r = t0; *v3_i = t1;

				mc_zunit1x2l(v0_r, v0_i, v2_r, v2_i);
				mc_zunit1x2l(v1_r, v1_i, v3_r, v3_i);
			} else {
				mc_zeye2x2l(
					  v0_r, v0_i, v1_r, v1_i
					, v2_r, v2_i, v3_r, v3_i
					, 0
				);
			}
		}
		return a[2] == a[1] ? 1 : (t2 > 0.0L ? 0 : 1);
	}
	return -1;
}

#endif /* !MC_ZREIG2X2_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_eigsy2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fmax.h>
#include <macadam/details/math/mc_fmin.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/math/mc_zsqrt.h>
#include <macadam/details/numa/mc_det2x2.h>
#include <macadam/details/numa/mc_trace2x2.h>

#ifndef MC_EIGSY2X2_H
#define MC_EIGSY2X2_H

#pragma mark - mc_eigsy2x2 -

MC_TARGET_FUNC int mc_eigsy2x2f(const float a[4], float e[2], float * v)
{
	const int wantv = mc_nonnullptr(v);
	float t0        = 0.0f, t1, t2;

	if (a[1] == a[2]) {
		if (mc_fabsf(a[2]) == 0.0f) {
			e[0] = mc_fminf(a[0], a[3]);
			e[1] = mc_fmaxf(a[0], a[3]);
			if (wantv) {
				mc_eye2x2f(v, 0);
			}
			return -3;
		} else {
			const float trc = mc_trace2x2f(a);
			const float det = mc_det2x2f(a);
			t1              = mc_raise2f(trc) - 4.0f * det;

			if (t1 >= 0.0f) {
				t1 = mc_sqrtf(t1);
			} else {
				mc_zsqrtf(&t0, &t1, t1, 0.0f);
			}
			t0 = -trc;
			t0 = -0.5f * (t0 + (t0 > 0.0f ? 1.0f : -1.0f) * t1);
			t1 = det / t0;

			if (mc_fabsf(t0) > mc_fabsf(t1)) {
				t2 = t0;
				t0 = t1;
				t1 = t2;
			}
			if (mc_fabsf(t0) < MCLIMITS_EPSILONF || mc_fabsf(t1) < MCLIMITS_EPSILONF) {
				e[0] = 1.0f; e[1] = 1.0f;
				if (wantv) {
					mc_eye2x2f(v, 0);
				}
				return -2;
			}
			e[0] = t0; e[1] = t1;
			if (wantv) {
				t0 =  a[3] - t0;
				t1 = -a[2];
				t2 = mc_sqrtf(mc_raise2f(t0) + mc_raise2f(t1));
				t0 = t0 / t2;
				t1 = t1 / t2;
				v[0] = t0; v[1] =  t1;
				v[2] = t1; v[3] = -t0;
			}
			return 0;
		}
	} else {
		e[0] = 1.0f; e[1] = 1.0f;
		if (wantv) {
			mc_eye2x2f(v, 0);
		}
	}
	return -1;
}

MC_TARGET_FUNC int mc_eigsy2x2ff(const float a[4], double e[2], double * v)
{
	const int wantv = mc_nonnullptr(v);
	double t0       = 0.0, t1, t2;

	if (a[1] == a[2]) {
		if (mc_fabsf(a[2]) == 0.0f) {
			e[0] = mc_cast(double, mc_fminf(a[0], a[3]));
			e[1] = mc_cast(double, mc_fmaxf(a[0], a[3]));
			if (wantv) {
				mc_eye2x2(v, 0);
			}
			return -3;
		} else {
			const double trc = mc_trace2x2ff(a);
			const double det = mc_det2x2ff(a);
			t1               = mc_raise2(trc) - 4.0 * det;

			if (t1 >= 0.0) {
				t1 = mc_sqrt(t1);
			} else {
				mc_zsqrt(&t0, &t1, t1, 0.0);
			}
			t0 = -trc;
			t0 = -0.5 * (t0 + (t0 > 0.0 ? 1.0 : -1.0) * t1);
			t1 = det / t0;

			if (mc_fabs(t0) > mc_fabs(t1)) {
				t2 = t0;
				t0 = t1;
				t1 = t2;
			}
			if (mc_fabs(t0) < MCLIMITS_EPSILON || mc_fabs(t1) < MCLIMITS_EPSILON) {
				e[0] = 1.0; e[1] = 1.0;
				if (wantv) {
					mc_eye2x2(v, 0);
				}
				return -2;
			}
			e[0] = t0; e[1] = t1;
			if (wantv) {
				t0 =  a[3] - t0;
				t1 = -a[2];
				t2 = mc_sqrt(mc_raise2(t0) + mc_raise2(t1));
				t0 = t0 / t2;
				t1 = t1 / t2;
				v[0] = t0; v[1] =  t1;
				v[2] = t1; v[3] = -t0;
			}
			return 0;
		}
	} else {
		e[0] = 1.0; e[1] = 1.0;
		if (wantv) {
			mc_eye2x2(v, 0);
		}
	}
	return -1;
}

MC_TARGET_FUNC int mc_eigsy2x2(const double a[4], double e[2], double * v)
{
	const int wantv = mc_nonnullptr(v);
	double t0       = 0.0, t1, t2;

	if (a[1] == a[2]) {
		if (mc_fabs(a[2]) == 0.0) {
			e[0] = mc_fmin(a[0], a[3]);
			e[1] = mc_fmax(a[0], a[3]);
			if (wantv) {
				mc_eye2x2(v, 0);
			}
			return -3;
		} else {
			const double trc = mc_trace2x2(a);
			const double det = mc_det2x2(a);
			t1               = mc_raise2(trc) - 4.0 * det;

			if (t1 >= 0.0) {
				t1 = mc_sqrt(t1);
			} else {
				mc_zsqrt(&t0, &t1, t1, 0.0);
			}
			t0 = -trc;
			t0 = -0.5 * (t0 + (t0 > 0.0 ? 1.0 : -1.0) * t1);
			t1 = det / t0;

			if (mc_fabs(t0) > mc_fabs(t1)) {
				t2 = t0;
				t0 = t1;
				t1 = t2;
			}
			if (mc_fabs(t0) < MCLIMITS_EPSILON || mc_fabs(t1) < MCLIMITS_EPSILON) {
				e[0] = 1.0; e[1] = 1.0;
				if (wantv) {
					mc_eye2x2(v, 0);
				}
				return -2;
			}
			e[0] = t0; e[1] = t1;
			if (wantv) {
				t0 =  a[3] - t0;
				t1 = -a[2];
				t2 = mc_sqrt(mc_raise2(t0) + mc_raise2(t1));
				t0 = t0 / t2;
				t1 = t1 / t2;
				v[0] = t0; v[1] =  t1;
				v[2] = t1; v[3] = -t0;
			}
			return 0;
		}
	} else {
		e[0] = 1.0; e[1] = 1.0;
		if (wantv) {
			mc_eye2x2(v, 0);
		}
	}
	return -1;
}

MC_TARGET_FUNC int mc_eigsy2x2l(const long double a[4], long double e[2], long double * v)
{
	const int wantv      = mc_nonnullptr(v);
	long double t0       = 0.0L, t1, t2;

	if (a[1] == a[2]) {
		if (mc_fabsl(a[2]) == 0.0L) {
			e[0] = mc_fminl(a[0], a[3]);
			e[1] = mc_fmaxl(a[0], a[3]);
			if (wantv) {
				mc_eye2x2l(v, 0);
			}
			return -3;
		} else {
			const long double trc = mc_trace2x2l(a);
			const long double det = mc_det2x2l(a);
			t1                    = mc_raise2l(trc) - 4.0L * det;

			if (t1 >= 0.0L) {
				t1 = mc_sqrtl(t1);
			} else {
				mc_zsqrtl(&t0, &t1, t1, 0.0L);
			}
			t0 = -trc;
			t0 = -0.5L * (t0 + (t0 > 0.0L ? 1.0L : -1.0L) * t1);
			t1 = det / t0;

			if (mc_fabsl(t0) > mc_fabsl(t1)) {
				t2 = t0;
				t0 = t1;
				t1 = t2;
			}
			if (mc_fabsl(t0) < MCLIMITS_EPSILONL || mc_fabsl(t1) < MCLIMITS_EPSILONL) {
				e[0] = 1.0L; e[1] = 1.0L;
				if (wantv) {
					mc_eye2x2l(v, 0);
				}
				return -2;
			}
			e[0] = t0; e[1] = t1;
			if (wantv) {
				t0   =  a[3] - t0;
				t1   = -a[2];
				t2   = mc_sqrtl(mc_raise2l(t0) + mc_raise2l(t1));
				t0   = t0 / t2;
				t1   = t1 / t2;
				v[0] = t0; v[1] =  t1;
				v[2] = t1; v[3] = -t0;
			}
			return 0;
		}
	} else {
		e[0] = 1.0L; e[1] = 1.0L;
		if (wantv) {
			mc_eye2x2l(v, 0);
		}
	}
	return -1;
}

#endif /* !MC_EIGSY2X2_H */

/* EOF */
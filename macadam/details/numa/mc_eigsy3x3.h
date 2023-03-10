//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_eigsy3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_hypot2.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/numa/mc_eye3x3.h>
#include <macadam/mcswap.h>

#ifndef MC_EIGSY3X3_H
#define MC_EIGSY3X3_H

#pragma mark - mc_tredsy3x3 -

MC_TARGET_PROC int mc_tredsy3x3f(const float a[9], float * q, float d[3], float e[2])
{
//!# Close-formish expression by rref.
	const float tol = 2.0f * MCLIMITS_EPSILONF;

	int wantq = mc_nonnullptr(q);
	float mag, s;

	float a11 = a[0], a12 = a[1], a13 = a[2];
	float             a22 = a[4], a23 = a[5];
	float                         a33 = a[8];

	if (wantq) {
		mc_eye3x3f(q, 0);
	}

	d[0] = a11; d[1] = a22; d[2] = a33;
	e[0] = a12; e[1] = a23;
	if (mc_fabsf(a13) >= tol && !(a12 == 0.0f && a13 == 0.0f)) {
		mag = mc_hypot2f(a12, a13);
		if (mc_fabsf(mag) >= tol) {
			if (a11 > 0.0f) {
				mag = -mag;
			}
			a12  = a12 / mag;
			a13  = a13 / mag;
			s    = (2.0f * a12 * a23) + a13 * (a33 - a22);
			d[1] = a22 + a13 * s; d[2] = a33 - a13 * s;
			e[0] = mag; e[1] = a23 - a12 * s;
			if (wantq) {
				q[4] = a12; q[5] = a13;
				q[7] = a13; q[8] =-a12;
			}
			return 0;
		}
	}
	return -1;
}

MC_TARGET_PROC int mc_tredsy3x3(const double a[9], double * q, double d[3], double e[2])
{
//!# Close-formish expression by rref.
	const double tol = MCLIMITS_TINY;

	int wantq = mc_nonnullptr(q);
	double mag, s;

	double a11 = a[0], a12 = a[1], a13 = a[2];
	double             a22 = a[4], a23 = a[5];
	double                         a33 = a[8];

	if (wantq) {
		mc_eye3x3(q, 0);
	}

	d[0] = a11; d[1] = a22; d[2] = a33;
	e[0] = a12; e[1] = a23;
	if (mc_fabs(a13) >= tol && !(a12 == 0.0 && a13 == 0.0)) {
		mag = mc_sqrt(mc_raise2(a12) + mc_raise2(a13));
		if (mc_fabs(mag) >= tol) {
			if (a11 > 0.0) {
				mag = -mag;
			}
			a12  = a12 / mag;
			a13  = a13 / mag;
			s    = 2.0 * a12 * a23 + a13 * (a33 - a22);
			d[1] = a22 + a13 * s; d[2] = a33 - a13 * s;
			e[0] = mag; e[1] = a23 - a12 * s;
			if (wantq) {
				q[4] = a12; q[5] = a13;
				q[7] = a13; q[8] =-a12;
			}
			return 0;
		}
	}
	return -1;
}

MC_TARGET_PROC int mc_tredsy3x3l(const long double a[9], long double * q, long double d[3], long double e[2])
{
//!# Close-formish expression by rref.
	const long double tol = 2.0L * MCLIMITS_EPSILONL;

	int wantq = mc_nonnullptr(q);
	long double mag, s;

	long double a11 = a[0], a12 = a[1], a13 = a[2];
	long double             a22 = a[4], a23 = a[5];
	long double                         a33 = a[8];

	if (wantq) {
		mc_eye3x3l(q, 0);
	}

	d[0] = a11; d[1] = a22; d[2] = a33;
	e[0] = a12; e[1] = a23;
	if (mc_fabsl(a13) >= tol && !(a12 == 0.0L && a13 == 0.0L)) {
		mag = mc_sqrtl(mc_raise2l(a12) + mc_raise2l(a13));
		if (mc_fabsl(mag) >= tol) {
			if (a11 > 0.0) {
				mag = -mag;
			}
			a12  = a12 / mag;
			a13  = a13 / mag;
			s    = 2.0L * a12 * a23 + a13 * (a33 - a22);
			d[1] = a22 + a13 * s; d[2] = a33 - a13 * s;
			e[0] = mag; e[1] = a23 - a12 * s;
			if (wantq) {
				q[4] = a12; q[5] = a13;
				q[7] = a13; q[8] =-a12;
			}
			return 0;
		}
	}
	return -1;
}

#pragma mark - mc_tredql3x3 -

MC_TARGET_PROC int mc_tredql3x3f(float * a, float d[3], float e[3])
{
	const int wanta = mc_nonnullptr(a);
	int z           = 0;
	float b, c, f, h, p, r, s, t;
	int i, j, k, n;

	for (k = 0; k < 2; k++) {
		n = 0;
		do {
			for (j = k; j <= 1; j++) {
				h = mc_fabsf(d[j]) + mc_fabsf(d[j + 1]);
				if (mc_fabsf(e[j]) + h == h) {
					break;
				}
			}

			if (j == k) {
				break;
			}
			z = !z;

			h = (d[k + 1] - d[k]) / (e[k] + e[k]);
			r = mc_hypot2f(1.0f, h);
			h = h > 0.0f ? (d[j] - d[k] + e[k] / (h + r)) : (d[j] - d[k] + e[k] / (h - r));

			c = 1.0f;
			s = 1.0f;
			p = 0.0f;

			for (i = (j - 1); i >= k; i--) {
				f = s * e[i];
				b = c * e[i];
				if (mc_fabsf(f) > mc_fabsf(h)) {
					c        = h / f;
					r        = mc_hypot2f(1.0f, c);
					e[i + 1] = f * r;
					s        = 1.0f / r;
					c        = c * s;
				} else {
					s        = f / h;
					r        = mc_hypot2f(1.0f, s);
					e[i + 1] = h * r;
					c        = 1.0f / r;
					s        = s * c;
				}

				h        = d[i + 1] - p;
				r        = (d[i] - h) * s + 2.0f * c * b;
				p        = s * r;
				d[i + 1] = h + p;
				h        = c * r - b;

				if (wanta) {
					t            = a[i + 1];
					a[i + 1    ] = s * a[i] + c * t;
					a[i        ] = c * a[i] - s * t;

					t            = a[3 + i + 1];
					a[3 + i + 1] = s * a[3 + i] + c * t;
					a[3 + i    ] = c * a[3 + i] - s * t;

					t            = a[6 + i + 1];
					a[6 + i + 1] = s * a[6 + i] + c * t;
					a[6 + i    ] = c * a[6 + i] - s * t;
				}
			}
			d[k] = d[k] - p;
			e[k] = h;
			e[j] = 0.0f;
		} while (++n <= 30);
		if (n >= 30) {
			return -1;
		}
	}

	if (mc_fabsf(d[0]) > mc_fabsf(d[1])) {
		mcswap_var(e[0], d[0], d[1]);
		if (wanta) {
			mcswap_var(e[0], a[0], a[1]);
			mcswap_var(e[0], a[3], a[4]);
			mcswap_var(e[0], a[6], a[7]);
		}
	}
	if (mc_fabsf(d[0]) > mc_fabsf(d[2])) {
		mcswap_var(e[0], d[0], d[2]);
		if (wanta) {
			mcswap_var(e[0], a[0], a[2]);
			mcswap_var(e[0], a[3], a[5]);
			mcswap_var(e[0], a[6], a[8]);
		}
	}
	if (mc_fabsf(d[1]) > mc_fabsf(d[2])) {
		mcswap_var(e[0], d[1], d[2]);
		if (wanta) {
			mcswap_var(e[0], a[1], a[2]);
			mcswap_var(e[0], a[4], a[5]);
			mcswap_var(e[0], a[7], a[8]);
		}
	}
	if (z && wanta) {
		a[2] = -a[2];
		a[5] = -a[5];
		a[8] = -a[8];
	}
	return 0;
}

MC_TARGET_PROC int mc_tredql3x3(double * a, double d[3], double e[3])
{
	const int wanta = mc_nonnullptr(a);
	int z           = 0;
	double b, c, f, h, p, r, s, t;
	int i, j, k, n;

	for (k = 0; k < 2; k++) {
		n = 0;
		do {
			for (j = k; j <= 1; j++) {
				h = mc_fabs(d[j]) + mc_fabs(d[j + 1]);
				if (mc_fabs(e[j]) + h == h) {
					break;
				}
			}

			if (j == k) {
				break;
			}
			z = !z;

			h = (d[k + 1] - d[k]) / (e[k] + e[k]);
			r = mc_hypot2(1.0, h);
			h = h > 0.0 ? (d[j] - d[k] + e[k] / (h + r)) : (d[j] - d[k] + e[k] / (h - r));

			c = 1.0;
			s = 1.0;
			p = 0.0;

			for (i = (j - 1); i >= k; i--) {
				f = s * e[i];
				b = c * e[i];
				if (mc_fabs(f) > mc_fabs(h)) {
					c        = h / f;
					r        = mc_hypot2(1.0, c);
					e[i + 1] = f * r;
					s        = 1.0 / r;
					c        = c * s;
				} else {
					s        = f / h;
					r        = mc_hypot2(1.0, s);
					e[i + 1] = h * r;
					c        = 1.0 / r;
					s        = s * c;
				}

				h        = d[i + 1] - p;
				r        = (d[i] - h) * s + 2.0 * c * b;
				p        = s * r;
				d[i + 1] = h + p;
				h        = c * r - b;

				if (wanta) {
					t            = a[i + 1];
					a[i + 1    ] = s * a[i] + c * t;
					a[i        ] = c * a[i] - s * t;

					t            = a[3 + i + 1];
					a[3 + i + 1] = s * a[3 + i] + c * t;
					a[3 + i    ] = c * a[3 + i] - s * t;

					t            = a[6 + i + 1];
					a[6 + i + 1] = s * a[6 + i] + c * t;
					a[6 + i    ] = c * a[6 + i] - s * t;
				}
			}
			d[k] = d[k] - p;
			e[k] = h;
			e[j] = 0.0;
		} while (++n <= 90);
		if (n >= 90) {
			return -1;
		}
	}

	if (mc_fabs(d[0]) > mc_fabs(d[1])) {
		mcswap_var(e[0], d[0], d[1]);
		if (wanta) {
			mcswap_var(e[0], a[0], a[1]);
			mcswap_var(e[0], a[3], a[4]);
			mcswap_var(e[0], a[6], a[7]);
		}
	}
	if (mc_fabs(d[0]) > mc_fabs(d[2])) {
		mcswap_var(e[0], d[0], d[2]);
		if (wanta) {
			mcswap_var(e[0], a[0], a[2]);
			mcswap_var(e[0], a[3], a[5]);
			mcswap_var(e[0], a[6], a[8]);
		}
	}
	if (mc_fabs(d[1]) > mc_fabs(d[2])) {
		mcswap_var(e[0], d[1], d[2]);
		if (wanta) {
			mcswap_var(e[0], a[1], a[2]);
			mcswap_var(e[0], a[4], a[5]);
			mcswap_var(e[0], a[7], a[8]);
		}
	}
	if (z && wanta) {
		a[2] = -a[2];
		a[5] = -a[5];
		a[8] = -a[8];
	}
	return 0;
}

MC_TARGET_PROC int mc_tredql3x3l(long double * a, long double d[3], long double e[3])
{
	const int wanta = mc_nonnullptr(a);
	int z           = 0;
	long double b, c, f, h, p, r, s, t;
	int i, j, k, n;

	for (k = 0; k < 2; k++) {
		n = 0;
		do {
			for (j = k; j <= 1; j++) {
				h = mc_fabsl(d[j]) + mc_fabsl(d[j + 1]);
				if (mc_fabsl(e[j]) + h == h) {
					break;
				}
			}

			if (j == k) {
				break;
			}
			z = !z;

			h = (d[k + 1] - d[k]) / (e[k] + e[k]);
			r = mc_hypot2l(1.0L, h);
			h = h > 0.0L ? (d[j] - d[k] + e[k] / (h + r)) : (d[j] - d[k] + e[k] / (h - r));

			c = 1.0L;
			s = 1.0L;
			p = 0.0L;

			for (i = (j - 1); i >= k; i--) {
				f = s * e[i];
				b = c * e[i];
				if (mc_fabsl(f) > mc_fabsl(h)) {
					c        = h / f;
					r        = mc_hypot2l(1.0L, c);
					e[i + 1] = f * r;
					s        = 1.0L / r;
					c        = c * s;
				} else {
					s        = f / h;
					r        = mc_hypot2l(1.0L, s);
					e[i + 1] = h * r;
					c        = 1.0L / r;
					s        = s * c;
				}

				h        = d[i + 1] - p;
				r        = (d[i] - h) * s + 2.0L * c * b;
				p        = s * r;
				d[i + 1] = h + p;
				h        = c * r - b;

				if (wanta) {
					t            = a[i + 1];
					a[i + 1    ] = s * a[i] + c * t;
					a[i        ] = c * a[i] - s * t;

					t            = a[3 + i + 1];
					a[3 + i + 1] = s * a[3 + i] + c * t;
					a[3 + i    ] = c * a[3 + i] - s * t;

					t            = a[6 + i + 1];
					a[6 + i + 1] = s * a[6 + i] + c * t;
					a[6 + i    ] = c * a[6 + i] - s * t;
				}
			}
			d[k] = d[k] - p;
			e[k] = h;
			e[j] = 0.0L;
		} while (++n <= 120);
		if (n >= 120) {
			return -1;
		}
	}

	if (mc_fabsl(d[0]) > mc_fabsl(d[1])) {
		mcswap_var(e[0], d[0], d[1]);
		if (wanta) {
			mcswap_var(e[0], a[0], a[1]);
			mcswap_var(e[0], a[3], a[4]);
			mcswap_var(e[0], a[6], a[7]);
		}
	}
	if (mc_fabsl(d[0]) > mc_fabsl(d[2])) {
		mcswap_var(e[0], d[0], d[2]);
		if (wanta) {
			mcswap_var(e[0], a[0], a[2]);
			mcswap_var(e[0], a[3], a[5]);
			mcswap_var(e[0], a[6], a[8]);
		}
	}
	if (mc_fabsl(d[1]) > mc_fabsl(d[2])) {
		mcswap_var(e[0], d[1], d[2]);
		if (wanta) {
			mcswap_var(e[0], a[1], a[2]);
			mcswap_var(e[0], a[4], a[5]);
			mcswap_var(e[0], a[7], a[8]);
		}
	}
	if (z && wanta) {
		a[2] = -a[2];
		a[5] = -a[5];
		a[8] = -a[8];
	}
	return 0;
}

#pragma mark - mc_eigsyq3x3 -

MC_TARGET_FUNC int mc_eigsyq3x3f(const float a[9], float e[3], float * v)
{
	int r;
	float w[3] = { 0 };
	if (0 == (r = mc_tredsy3x3f(a, v, e, w))) {
		r = mc_tredql3x3f(v, e, w);
	}
	return r;
}

MC_TARGET_FUNC int mc_eigsyq3x3(const double a[9], double e[3], double * v)
{
	int r;
	double w[3] = { 0 };
	if (0 == (r = mc_tredsy3x3(a, v, e, w))) {
		r = mc_tredql3x3(v, e, w);
	}
	return r;
}

MC_TARGET_FUNC int mc_eigsyq3x3l(const long double a[9], long double e[3], long double * v)
{
	int r;
	long double w[3] = { 0 };
	if (0 == (r = mc_tredsy3x3l(a, v, e, w))) {
		r = mc_tredql3x3l(v, e, w);
	}
	return r;
}

#pragma mark - mc_eigsy3x3 -

MC_TARGET_FUNC int mc_eigsy3x3f(const float a[9], float e[3], float * v)
{
	const int wantv = mc_nonnullptr(v);
//!# Number of Jacobi iterations.
	int i           = 0;
//!# Too low values guard.
	const float m   = MCLIMITS_TINYF;
//!# Max number of iteration for convergence.
	const int j     = 30;

//!# Copying upper triangle of the given symmetric system.
	float a11 = a[0], a12 = a[1], a13 = a[2];
	float             a22 = a[4], a23 = a[5];
	float                         a33 = a[8];

//!# Absolute values of off-diagonal items.
	float aa12 = mc_fabsf(a12);
	float aa13 = mc_fabsf(a13);
	float aa23 = mc_fabsf(a23);

//!# Working vectors. 
	float v11 = 1.0f, v12 = 0.0f, v13 = 0.0f;
	float v21 = 0.0f, v22 = 1.0f, v23 = 0.0f;
	float v31 = 0.0f, v32 = 0.0f, v33 = 1.0f;

//!# Jacobi rotation variables.
	float c, r, s, t, u;
	float apr, aqr, vpr, vqr;

//!# Applying Jacobi rotations until all off-diagonal items are 0.
	for (; aa12 + aa13 + aa23 > 0.0f; i++) {
		if (i >= j) {
			e[0] = 1.0f, e[1] = 1.0f, e[2] = 1.0f;
			if (wantv) {
				mc_eye3x3f(v, 0);
			}
			break;
		}
		if (aa12 >= aa13 && aa12 >= aa23) {
			u = a22 - a11;
			if (mc_fabsf(a12) < m * mc_fabsf(u)) {
				t = a12 / u;
			} else {
				r = 0.5f * u / a12;
				t = ((r >= 0.0f) ? 1.0f / (r + mc_hypot2f(1.0f, r)) : 1.0f / (r - mc_hypot2f(1.0f, r)));
			}
			c   = 1.0f / mc_hypot2f(1.0f, t);
			s   = t * c;
			u   = s / (1.0f + c);
			r   = t * a12;
			a11 = a11 - r;
			a22 = a22 + r;
			a12 = 0.0f;
			apr = a13;
			aqr = a23;
			a13 = apr - s * (aqr + apr * u);
			a23 = aqr + s * (apr - aqr * u);
			vpr = v11;
			vqr = v12;
			v11 = vpr - s * (vqr + vpr * u);
			v12 = vqr + s * (vpr - vqr * u);
			vpr = v21;
			vqr = v22;
			v21 = vpr - s * (vqr + vpr * u);
			v22 = vqr + s * (vpr - vqr * u);
			vpr = v31;
			vqr = v32;
			v31 = vpr - s * (vqr + vpr * u);
			v32 = vqr + s * (vpr - vqr * u);
		} else if (aa13 >= aa12 && aa13 >= aa23) {
			u = a33 - a11;
			if (mc_fabsf(a13) < m * mc_fabsf(u)) {
				t = a13 / u;
			} else {
				r = 0.5f * u / a13;
				t = ((r >= 0.0f) ? 1.0f / (r + mc_hypot2f(1.0f, r)) : 1.0f / (r - mc_hypot2f(1.0f, r)));
			}
			c   = 1.0f / mc_hypot2f(1.0f, t);
			s   = t * c;
			u   = s / (1.0f + c);
			r   = t * a13;
			a11 = a11 - r;
			a33 = a33 + r;
			a13 = 0.0f;
			apr = a12;
			aqr = a23;
			a12 = apr - s * (aqr + apr * u);
			a23 = aqr + s * (apr - aqr * u);
			vpr = v11;
			vqr = v13;
			v11 = vpr - s * (vqr + vpr * u);
			v13 = vqr + s * (vpr - vqr * u);
			vpr = v21;
			vqr = v23;
			v21 = vpr - s * (vqr + vpr * u);
			v23 = vqr + s * (vpr - vqr * u);
			vpr = v31;
			vqr = v33;
			v31 = vpr - s * (vqr + vpr * u);
			v33 = vqr + s * (vpr - vqr * u);
		} else {
			u = a33 - a22;
			if (mc_fabsf(a23) < m * mc_fabsf(u)) {
				t = a23 / u;
			} else {
				r = 0.5f * u / a23;
				t = ((r >= 0.0f) ? 1.0f / (r + mc_hypot2f(1.0f, r)) : 1.0f / (r - mc_hypot2f(1.0f, r)));
			}
			c   = 1.0f / mc_hypot2f(1.0f, t);
			s   = t * c;
			u   = s / (1.0f + c);
			r   = t * a23;
			a22 = a22 - r;
			a33 = a33 + r;
			a23 = 0.0f;
			apr = a12;
			aqr = a13;
			a12 = apr - s * (aqr + apr * u);
			a13 = aqr + s * (apr - aqr * u);
			vpr = v12;
			vqr = v13;
			v12 = vpr - s * (vqr + vpr * u);
			v13 = vqr + s * (vpr - vqr * u);
			vpr = v22;
			vqr = v23;
			v22 = vpr - s * (vqr + vpr * u);
			v23 = vqr + s * (vpr - vqr * u);
			vpr = v32;
			vqr = v33;
			v32 = vpr - s * (vqr + vpr * u);
			v33 = vqr + s * (vpr - vqr * u);
		}
		if (mc_fabsf(v11) < MCLIMITS_EPSILONF) {
			v11 = 0.0f;
		}
		if (mc_fabsf(v12) < MCLIMITS_EPSILONF) {
			v12 = 0.0f;
		}
		if (mc_fabsf(v13) < MCLIMITS_EPSILONF) {
			v13 = 0.0f;
		}
		if (mc_fabsf(v21) < MCLIMITS_EPSILONF) {
			v21 = 0.0f;
		}
		if (mc_fabsf(v22) < MCLIMITS_EPSILONF) {
			v22 = 0.0f;
		}
		if (mc_fabsf(v23) < MCLIMITS_EPSILONF) {
			v23 = 0.0f;
		}
		if (mc_fabsf(v31) < MCLIMITS_EPSILONF) {
			v31 = 0.0f;
		}
		if (mc_fabsf(v32) < MCLIMITS_EPSILONF) {
			v32 = 0.0f;
		}
		if (mc_fabsf(v33) < MCLIMITS_EPSILONF) {
			v33 = 0.0f;
		}
		aa12 = mc_fabsf(a12);
		aa13 = mc_fabsf(a13);
		aa23 = mc_fabsf(a23);
		if (aa12 < MCLIMITS_EPSILONF) {
			aa12 = 0.0f;
		}
		if (aa13 < MCLIMITS_EPSILONF) {
			aa13 = 0.0f;
		}
		if (aa23 < MCLIMITS_EPSILONF) {
			aa23 = 0.0f;
		}
	}
	if (i < j) {
//!# Reordering eigenvalues and eigenvectors (absolute ascending i.e smaller first).
		if (mc_fabsf(a11) > mc_fabsf(a22)) {
			mcswap_var(t, a11, a22);
			if (wantv) {
				mcswap_var(t, v11, v12);
				mcswap_var(t, v21, v22);
				mcswap_var(t, v31, v32);
			}
		}
		if (mc_fabsf(a11) > mc_fabsf(a33)) {
			mcswap_var(t, a11, a33);
			if (wantv) {
				mcswap_var(t, v11, v13);
				mcswap_var(t, v21, v23);
				mcswap_var(t, v31, v33);
			}
		}
		if (mc_fabsf(a22) > mc_fabsf(a33)) {
			mcswap_var(t, a22, a33);
			if (wantv) {
				mcswap_var(t, v12, v13);
				mcswap_var(t, v22, v23);
				mcswap_var(t, v32, v33);
			}
		}
		e[0] = a11; e[1] = a22; e[2] = a33;
		if (wantv) {
			v[0] = v11; v[1] = v12; v[2] = v13;
			v[3] = v21; v[4] = v22; v[5] = v23;
			v[6] = v31; v[7] = v32; v[8] = v33;
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_eigsy3x3(const double a[9], double e[3], double * v)
{
	const int wantv = mc_nonnullptr(v);
//!# Number of Jacobi iterations.
	int i           = 0;
//!# Too low values guard.
	const double m  = MCLIMITS_TINY;
//!# Max number of iteration for convergence.
	const int j     = 90;

//!# Copying upper triangle of the given symmetric system.
	double a11 = a[0], a12 = a[1], a13 = a[2];
	double             a22 = a[4], a23 = a[5];
	double                         a33 = a[8];

//!# Absolute values of off-diagonal items.
	double aa12 = mc_fabs(a12);
	double aa13 = mc_fabs(a13);
	double aa23 = mc_fabs(a23);

//!# Working vectors. 
	double v11 = 1.0, v12 = 0.0, v13 = 0.0;
	double v21 = 0.0, v22 = 1.0, v23 = 0.0;
	double v31 = 0.0, v32 = 0.0, v33 = 1.0;

//!# Jacobi rotation variables.
	double c, r, s, t, u;
	double apr, aqr, vpr, vqr;

//!# Applying Jacobi rotations until all off-diagonal items are 0.
	for (; aa12 + aa13 + aa23 > 0.0; i++) {
		if (i >= j) {
			e[0] = 1.0, e[1] = 1.0, e[2] = 1.0;
			if (wantv) {
				mc_eye3x3(v, 0);
			}
			break;
		}
		if (aa12 >= aa13 && aa12 >= aa23) {
			u = a22 - a11;
			if (mc_fabs(a12) < m * mc_fabs(u)) {
				t = a12 / u;
			} else {
				r = 0.5 * u / a12;
				t = ((r >= 0.0) ? 1.0 / (r + mc_hypot2(1.0, r)) : 1.0 / (r - mc_hypot2(1.0, r)));
			}
			c   = 1.0 / mc_hypot2(1.0, t);
			s   = t * c;
			u   = s / (1.0 + c);
			r   = t * a12;
			a11 = a11 - r;
			a22 = a22 + r;
			a12 = 0.0;
			apr = a13;
			aqr = a23;
			a13 = apr - s * (aqr + apr * u);
			a23 = aqr + s * (apr - aqr * u);
			vpr = v11;
			vqr = v12;
			v11 = vpr - s * (vqr + vpr * u);
			v12 = vqr + s * (vpr - vqr * u);
			vpr = v21;
			vqr = v22;
			v21 = vpr - s * (vqr + vpr * u);
			v22 = vqr + s * (vpr - vqr * u);
			vpr = v31;
			vqr = v32;
			v31 = vpr - s * (vqr + vpr * u);
			v32 = vqr + s * (vpr - vqr * u);
		} else if (aa13 >= aa12 && aa13 >= aa23) {
			u = a33 - a11;
			if (mc_fabs(a13) < m * mc_fabs(u)) {
				t = a13 / u;
			} else {
				r = 0.5 * u / a13;
				t = ((r >= 0.0) ? 1.0 / (r + mc_hypot2(1.0, r)) : 1.0 / (r - mc_hypot2(1.0, r)));
			}
			c   = 1.0 / mc_hypot2(1.0, t);
			s   = t * c;
			u   = s / (1.0 + c);
			r   = t * a13;
			a11 = a11 - r;
			a33 = a33 + r;
			a13 = 0.0;
			apr = a12;
			aqr = a23;
			a12 = apr - s * (aqr + apr * u);
			a23 = aqr + s * (apr - aqr * u);
			vpr = v11;
			vqr = v13;
			v11 = vpr - s * (vqr + vpr * u);
			v13 = vqr + s * (vpr - vqr * u);
			vpr = v21;
			vqr = v23;
			v21 = vpr - s * (vqr + vpr * u);
			v23 = vqr + s * (vpr - vqr * u);
			vpr = v31;
			vqr = v33;
			v31 = vpr - s * (vqr + vpr * u);
			v33 = vqr + s * (vpr - vqr * u);
		} else {
			u = a33 - a22;
			if (mc_fabs(a23) < m * mc_fabs(u)) {
				t = a23 / u;
			} else {
				r = 0.5 * u / a23;
				t = ((r >= 0.0) ? 1.0 / (r + mc_hypot2(1.0, r)) : 1.0 / (r - mc_hypot2(1.0, r)));
			}
			c   = 1.0 / mc_hypot2(1.0, t);
			s   = t * c;
			u   = s / (1.0 + c);
			r   = t * a23;
			a22 = a22 - r;
			a33 = a33 + r;
			a23 = 0.0;
			apr = a12;
			aqr = a13;
			a12 = apr - s * (aqr + apr * u);
			a13 = aqr + s * (apr - aqr * u);
			vpr = v12;
			vqr = v13;
			v12 = vpr - s * (vqr + vpr * u);
			v13 = vqr + s * (vpr - vqr * u);
			vpr = v22;
			vqr = v23;
			v22 = vpr - s * (vqr + vpr * u);
			v23 = vqr + s * (vpr - vqr * u);
			vpr = v32;
			vqr = v33;
			v32 = vpr - s * (vqr + vpr * u);
			v33 = vqr + s * (vpr - vqr * u);
		}
		if (mc_fabs(v11) < MCLIMITS_EPSILON) {
			v11 = 0.0;
		}
		if (mc_fabs(v12) < MCLIMITS_EPSILON) {
			v12 = 0.0;
		}
		if (mc_fabs(v13) < MCLIMITS_EPSILON) {
			v13 = 0.0;
		}
		if (mc_fabs(v21) < MCLIMITS_EPSILON) {
			v21 = 0.0;
		}
		if (mc_fabs(v22) < MCLIMITS_EPSILON) {
			v22 = 0.0;
		}
		if (mc_fabs(v23) < MCLIMITS_EPSILON) {
			v23 = 0.0;
		}
		if (mc_fabs(v31) < MCLIMITS_EPSILON) {
			v31 = 0.0;
		}
		if (mc_fabs(v32) < MCLIMITS_EPSILON) {
			v32 = 0.0;
		}
		if (mc_fabs(v33) < MCLIMITS_EPSILON) {
			v33 = 0.0;
		}
		aa12 = mc_fabs(a12);
		aa13 = mc_fabs(a13);
		aa23 = mc_fabs(a23);
		if (aa12 < MCLIMITS_EPSILON) {
			aa12 = 0.0;
		}
		if (aa13 < MCLIMITS_EPSILON) {
			aa13 = 0.0;
		}
		if (aa23 < MCLIMITS_EPSILON) {
			aa23 = 0.0;
		}
	}
	if (i < j) {
//!# Reordering eigenvalues and eigenvectors (absolute ascending i.e smaller first).
		if (mc_fabs(a11) > mc_fabs(a22)) {
			mcswap_var(t, a11, a22);
			if (wantv) {
				mcswap_var(t, v11, v12);
				mcswap_var(t, v21, v22);
				mcswap_var(t, v31, v32);
			}
		}
		if (mc_fabs(a11) > mc_fabs(a33)) {
			mcswap_var(t, a11, a33);
			if (wantv) {
				mcswap_var(t, v11, v13);
				mcswap_var(t, v21, v23);
				mcswap_var(t, v31, v33);
			}
		}
		if (mc_fabs(a22) > mc_fabs(a33)) {
			mcswap_var(t, a22, a33);
			if (wantv) {
				mcswap_var(t, v12, v13);
				mcswap_var(t, v22, v23);
				mcswap_var(t, v32, v33);
			}
		}
		e[0] = a11; e[1] = a22; e[2] = a33;
		if (wantv) {
			v[0] = v11; v[1] = v12; v[2] = v13;
			v[3] = v21; v[4] = v22; v[5] = v23;
			v[6] = v31; v[7] = v32; v[8] = v33;
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_eigsy3x3l(const long double a[9], long double e[3], long double * v)
{
	const int wantv     = mc_nonnullptr(v);
//!# Number of Jacobi iterations.
	int i               = 0;
//!# Too low values guard.
	const long double m = MCLIMITS_TINYL;
//!# Max number of iteration for convergence.
	const int j         = 120;

//!# Copying upper triangle of the given symmetric system.
	long double a11 = a[0], a12 = a[1], a13 = a[2];
	long double             a22 = a[4], a23 = a[5];
	long double                         a33 = a[8];

//!# Absolute values of off-diagonal items.
	long double aa12 = mc_fabsl(a12);
	long double aa13 = mc_fabsl(a13);
	long double aa23 = mc_fabsl(a23);

//!# Working vectors. 
	long double v11 = 1.0L, v12 = 0.0L, v13 = 0.0L;
	long double v21 = 0.0L, v22 = 1.0L, v23 = 0.0L;
	long double v31 = 0.0L, v32 = 0.0L, v33 = 1.0L;

//!# Jacobi rotation variables.
	long double c, r, s, t, u;
	long double apr, aqr, vpr, vqr;

//!# Applying Jacobi rotations until all off-diagonal items are 0.
	for (; aa12 + aa13 + aa23 > 0.0L; i++) {
		if (i >= j) {
			e[0] = 1.0L, e[1] = 1.0L, e[2] = 1.0L;
			if (wantv) {
				mc_eye3x3l(v, 0);
			}
			break;
		}
		if (aa12 >= aa13 && aa12 >= aa23) {
			u = a22 - a11;
			if (mc_fabsl(a12) < m * mc_fabsl(u)) {
				t = a12 / u;
			} else {
				r = 0.5L * u / a12;
				t = ((r >= 0.0L) ? 1.0L / (r + mc_hypot2l(1.0L, r)) : 1.0L / (r - mc_hypot2l(1.0L, r)));
			}
			c   = 1.0L / mc_hypot2l(1.0L, t);
			s   = t * c;
			u   = s / (1.0L + c);
			r   = t * a12;
			a11 = a11 - r;
			a22 = a22 + r;
			a12 = 0.0L;
			apr = a13;
			aqr = a23;
			a13 = apr - s * (aqr + apr * u);
			a23 = aqr + s * (apr - aqr * u);
			vpr = v11;
			vqr = v12;
			v11 = vpr - s * (vqr + vpr * u);
			v12 = vqr + s * (vpr - vqr * u);
			vpr = v21;
			vqr = v22;
			v21 = vpr - s * (vqr + vpr * u);
			v22 = vqr + s * (vpr - vqr * u);
			vpr = v31;
			vqr = v32;
			v31 = vpr - s * (vqr + vpr * u);
			v32 = vqr + s * (vpr - vqr * u);
		} else if (aa13 >= aa12 && aa13 >= aa23) {
			u = a33 - a11;
			if (mc_fabsl(a13) < m * mc_fabsl(u)) {
				t = a13 / u;
			} else {
				r = 0.5L * u / a13;
				t = ((r >= 0.0L) ? 1.0L / (r + mc_hypot2l(1.0L, r)) : 1.0L / (r - mc_hypot2l(1.0L, r)));
			}
			c   = 1.0L / mc_hypot2l(1.0L, t);
			s   = t * c;
			u   = s / (1.0L + c);
			r   = t * a13;
			a11 = a11 - r;
			a33 = a33 + r;
			a13 = 0.0L;
			apr = a12;
			aqr = a23;
			a12 = apr - s * (aqr + apr * u);
			a23 = aqr + s * (apr - aqr * u);
			vpr = v11;
			vqr = v13;
			v11 = vpr - s * (vqr + vpr * u);
			v13 = vqr + s * (vpr - vqr * u);
			vpr = v21;
			vqr = v23;
			v21 = vpr - s * (vqr + vpr * u);
			v23 = vqr + s * (vpr - vqr * u);
			vpr = v31;
			vqr = v33;
			v31 = vpr - s * (vqr + vpr * u);
			v33 = vqr + s * (vpr - vqr * u);
		} else {
			u = a33 - a22;
			if (mc_fabsl(a23) < m * mc_fabsl(u)) {
				t = a23 / u;
			} else {
				r = 0.5L * u / a23;
				t = ((r >= 0.0L) ? 1.0L / (r + mc_hypot2l(1.0L, r)) : 1.0L / (r - mc_hypot2l(1.0L, r)));
			}
			c   = 1.0L / mc_hypot2l(1.0L, t);
			s   = t * c;
			u   = s / (1.0L + c);
			r   = t * a23;
			a22 = a22 - r;
			a33 = a33 + r;
			a23 = 0.0L;
			apr = a12;
			aqr = a13;
			a12 = apr - s * (aqr + apr * u);
			a13 = aqr + s * (apr - aqr * u);
			vpr = v12;
			vqr = v13;
			v12 = vpr - s * (vqr + vpr * u);
			v13 = vqr + s * (vpr - vqr * u);
			vpr = v22;
			vqr = v23;
			v22 = vpr - s * (vqr + vpr * u);
			v23 = vqr + s * (vpr - vqr * u);
			vpr = v32;
			vqr = v33;
			v32 = vpr - s * (vqr + vpr * u);
			v33 = vqr + s * (vpr - vqr * u);
		}
		if (mc_fabsl(v11) < MCLIMITS_EPSILONL) {
			v11 = 0.0L;
		}
		if (mc_fabsl(v12) < MCLIMITS_EPSILONL) {
			v12 = 0.0L;
		}
		if (mc_fabsl(v13) < MCLIMITS_EPSILONL) {
			v13 = 0.0L;
		}
		if (mc_fabsl(v21) < MCLIMITS_EPSILONL) {
			v21 = 0.0L;
		}
		if (mc_fabsl(v22) < MCLIMITS_EPSILONL) {
			v22 = 0.0L;
		}
		if (mc_fabsl(v23) < MCLIMITS_EPSILONL) {
			v23 = 0.0L;
		}
		if (mc_fabsl(v31) < MCLIMITS_EPSILONL) {
			v31 = 0.0L;
		}
		if (mc_fabsl(v32) < MCLIMITS_EPSILONL) {
			v32 = 0.0L;
		}
		if (mc_fabsl(v33) < MCLIMITS_EPSILONL) {
			v33 = 0.0L;
		}
		aa12 = mc_fabsl(a12);
		aa13 = mc_fabsl(a13);
		aa23 = mc_fabsl(a23);
		if (aa12 < MCLIMITS_EPSILONL) {
			aa12 = 0.0L;
		}
		if (aa13 < MCLIMITS_EPSILONL) {
			aa13 = 0.0L;
		}
		if (aa23 < MCLIMITS_EPSILONL) {
			aa23 = 0.0L;
		}
	}
	if (i < j) {
//!# Reordering eigenvalues and eigenvectors (absolute ascending i.e smaller first).
		if (mc_fabsl(a11) > mc_fabsl(a22)) {
			mcswap_var(t, a11, a22);
			if (wantv) {
				mcswap_var(t, v11, v12);
				mcswap_var(t, v21, v22);
				mcswap_var(t, v31, v32);
			}
		}
		if (mc_fabsl(a11) > mc_fabsl(a33)) {
			mcswap_var(t, a11, a33);
			if (wantv) {
				mcswap_var(t, v11, v13);
				mcswap_var(t, v21, v23);
				mcswap_var(t, v31, v33);
			}
		}
		if (mc_fabsl(a22) > mc_fabsl(a33)) {
			mcswap_var(t, a22, a33);
			if (wantv) {
				mcswap_var(t, v12, v13);
				mcswap_var(t, v22, v23);
				mcswap_var(t, v32, v33);
			}
		}
		e[0] = a11; e[1] = a22; e[2] = a33;
		if (wantv) {
			v[0] = v11; v[1] = v12; v[2] = v13;
			v[3] = v21; v[4] = v22; v[5] = v23;
			v[6] = v31; v[7] = v32; v[8] = v33;
		}
		return 0;
	}
	return -1;
}

#endif /* !MC_EIGSY3X3_H */

/* EOF */
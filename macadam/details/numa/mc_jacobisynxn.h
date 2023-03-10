//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_jacobisynxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_hypot2.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/numa/mc_copynxn.h>
#include <macadam/details/numa/mc_eyenxn.h>
#include <macadam/details/numa/mc_triussqrnxn.h>

#ifndef MC_JACOBISYNXN_H
#define MC_JACOBISYNXN_H

#pragma mark - mc_jacobisynxn -

MC_TARGET_FUNC int mc_jacobisynxnf(const int n, const float * a, float tol, float * e, float * v)
{
//!# Requires a[n x n], e[n x n] and v[n x n] if !null where 1 < n.
//!# A and E may be the same. Computing eigenvalues and right-hand
//!# side eigenvectors of a real symmetric matrix by Jacobi method.
	int i, j, k, l = 0;
	float c, s, scale, sumsq, diff, w, t;

	const int max   = 90;
	const int wantv = mc_nonnullptr(v);

	if (a != e) {
		mc_copynxnf(n, e, a);
	}
	if (wantv) {
		mc_eyenxnf(n, v, 0);
	}

	if (tol <= 0.0f) {
		tol = mc_fabsf(tol);
	}
	if (tol == 0.0f) {
		tol = MCLIMITS_TINYF;
	}

	while (1) {
//!# Computing the partial Frobenius norm, iterating the
//!# upper-triangle excluding the main-diagonal elements.
		mc_triussqrnxnf(n, e, &sumsq, &scale, 1);
		if (mc_fisnearf(scale * mc_sqrtf(sumsq * 2.0f), 0.0f, 1)) {
			break;
		}
		for (i = 0; i < n; i++) {
			for (j = i + 1; j < n; j++) {
				t = e[(n * i) + j];
				c = 1.0f;
				s = 0.0f;
//!# Forming symmetric Schur.
				if (mc_fabsf(t) > tol) {
					diff = (e[(n * j) + j] - e[(n * i) + i]) / (2.0f * t);
					t    = (mc_copysignf(1.0f, diff) > 0.0f
						? ( 1.0f / ( diff + mc_hypot2f(1.0f, diff)))
						: (-1.0f / (-diff + mc_hypot2f(1.0f, diff)))
					);
					w = 1.0f / mc_hypot2f(1.0f, t);
					c = w;
					s = t * w;
				}
//!# Rotating left.
				for (k = 0; k < n; k++) {
					t              = e[(n * i) + k];
					w              = e[(n * j) + k];
					e[(n * i) + k] = t * c - w * s;
					e[(n * j) + k] = t * s + w * c;
				}
//!# Rotating right.
				for (k = 0; k < n; k++) {
					t              = e[(n * k) + i];
					w              = e[(n * k) + j];
					e[(n * k) + i] = t * c - w * s;
					e[(n * k) + j] = t * s + w * c;
					if (wantv) {
						t              = v[(n * k) + i];
						w              = v[(n * k) + j];
						v[(n * k) + i] = t * c - w * s;
						v[(n * k) + j] = t * s + w * c;
					}
				}
			}
		}
		++l;
		if (l >= max) {
			return -1;
		}
	}

//!# Reordering eigenvalues and eigenvectors (absolute ascending i.e smaller first).
	for (i = 0; i < (n - 1); i++) {
		k = i;
		w = e[(n * k) + k];
		for (j = (i + 1); j < n; j++) {
			if (mc_fabsf(e[(n * j) + j]) < mc_fabsf(w)) {
				k = j;
				w = e[(n * k) + k];
			}
		}
		if (k > i) {
			e[(n * k) + k] = e[(n * i) + i];
			e[(n * i) + i] = w;
			if (wantv) {
				for (j = 0; j < n; j++) {
					w              = v[(n * j) + k];
					v[(n * j) + k] = v[(n * j) + i];
					v[(n * j) + i] = w;
				}
			}
		}
	}

	return 0;
}

MC_TARGET_FUNC int mc_jacobisynxnff(const int n, const float * a, float tol, double * e, double * v)
{
//!# Requires a[n x n], e[n x n] and v[n x n] if !null where 1 < n.
//!# Computing eigenvalues and right-hand side eigenvectors of a
//!# real symmetric matrix by Jacobi method.
	int i, j, k, l = 0;
	double c, s, scale, sumsq, diff, w, t, told;

	const int max   = 120;
	const int wantv = mc_nonnullptr(v);

	mc_copynxnff(n, e, a);

	if (wantv) {
		mc_eyenxn(n, v, 0);
	}

	if (tol <= 0.0f) {
		tol = mc_fabsf(tol);
	}
	if (tol == 0.0f) {
		tol = MCLIMITS_TINYF;
	}
	told = mc_cast(double, tol);

	while (1) {
//!# Computing the partial Frobenius norm, iterating the
//!# upper-triangle excluding the main-diagonal elements.
		mc_triussqrnxn(n, e, &sumsq, &scale, 1);
		if (mc_fisnear(scale * mc_sqrt(sumsq * 2.0), 0.0, 1)) {
			break;
		}
		for (i = 0; i < n; i++) {
			for (j = i + 1; j < n; j++) {
				t = e[(n * i) + j];
				c = 1.0;
				s = 0.0;
//!# Forming symmetric Schur.
				if (mc_fabs(t) > told) {
					diff = (e[(n * j) + j] - e[(n * i) + i]) / (2.0 * t);
					t    = (mc_copysign(1.0, diff) > 0.0
						? ( 1.0 / ( diff + mc_hypot2(1.0, diff)))
						: (-1.0 / (-diff + mc_hypot2(1.0, diff)))
					);
					w = 1.0 / mc_hypot2(1.0, t);
					c = w;
					s = t * w;
				}
//!# Rotating left.
				for (k = 0; k < n; k++) {
					t              = e[(n * i) + k];
					w              = e[(n * j) + k];
					e[(n * i) + k] = t * c - w * s;
					e[(n * j) + k] = t * s + w * c;
				}
//!# Rotating right.
				for (k = 0; k < n; k++) {
					t              = e[(n * k) + i];
					w              = e[(n * k) + j];
					e[(n * k) + i] = t * c - w * s;
					e[(n * k) + j] = t * s + w * c;
					if (wantv) {
						t              = v[(n * k) + i];
						w              = v[(n * k) + j];
						v[(n * k) + i] = t * c - w * s;
						v[(n * k) + j] = t * s + w * c;
					}
				}
			}
		}
		++l;
		if (l >= max) {
			return -1;
		}
	}

//!# Reordering eigenvalues and eigenvectors (absolute ascending i.e smaller first).
	for (i = 0; i < (n - 1); i++) {
		k = i;
		w = e[(n * k) + k];
		for (j = (i + 1); j < n; j++) {
			if (mc_fabs(e[(n * j) + j]) < mc_fabs(w)) {
				k = j;
				w = e[(n * k) + k];
			}
		}
		if (k > i) {
			e[(n * k) + k] = e[(n * i) + i];
			e[(n * i) + i] = w;
			if (wantv) {
				for (j = 0; j < n; j++) {
					w              = v[(n * j) + k];
					v[(n * j) + k] = v[(n * j) + i];
					v[(n * j) + i] = w;
				}
			}
		}
	}

	return 0;
}

MC_TARGET_FUNC int mc_jacobisynxn(const int n, const double * a, double tol, double * e, double * v)
{
//!# Requires a[n x n], e[n x n] and v[n x n] if !null where 1 < n.
//!# A and E may be the same. Computing eigenvalues and right-hand
//!# side eigenvectors of a real symmetric matrix by Jacobi method.
	int i, j, k, l = 0;
	double c, s, scale, sumsq, diff, w, t;

	const int max   = 120;
	const int wantv = mc_nonnullptr(v);

	if (a != e) {
		mc_copynxn(n, e, a);
	}

	if (wantv) {
		mc_eyenxn(n, v, 0);
	}

	if (tol <= 0.0) {
		tol = mc_fabs(tol);
	}
	if (tol == 0.0) {
		tol = MCLIMITS_TINY;
	}

	while (1) {
//!# Computing the partial Frobenius norm, iterating the
//!# upper-triangle excluding the main-diagonal elements.
		mc_triussqrnxn(n, e, &sumsq, &scale, 1);
		if (mc_fisnear(scale * mc_sqrt(sumsq * 2.0), 0.0, 1)) {
			break;
		}
		for (i = 0; i < n; i++) {
			for (j = i + 1; j < n; j++) {
				t = e[(n * i) + j];
				c = 1.0;
				s = 0.0;
//!# Forming symmetric Schur.
				if (mc_fabs(t) > tol) {
					diff = (e[(n * j) + j] - e[(n * i) + i]) / (2.0 * t);
					t    = (mc_copysign(1.0, diff) > 0.0
						? ( 1.0 / ( diff + mc_hypot2(1.0, diff)))
						: (-1.0 / (-diff + mc_hypot2(1.0, diff)))
					);
					w = 1.0 / mc_hypot2(1.0, t);
					c = w;
					s = t * w;
				}
//!# Rotating left.
				for (k = 0; k < n; k++) {
					t              = e[(n * i) + k];
					w              = e[(n * j) + k];
					e[(n * i) + k] = t * c - w * s;
					e[(n * j) + k] = t * s + w * c;
				}
//!# Rotating right.
				for (k = 0; k < n; k++) {
					t              = e[(n * k) + i];
					w              = e[(n * k) + j];
					e[(n * k) + i] = t * c - w * s;
					e[(n * k) + j] = t * s + w * c;
					if (wantv) {
						t              = v[(n * k) + i];
						w              = v[(n * k) + j];
						v[(n * k) + i] = t * c - w * s;
						v[(n * k) + j] = t * s + w * c;
					}
				}
			}
		}
		++l;
		if (l >= max) {
			return -1;
		}
	}

//!# Reordering eigenvalues and eigenvectors (absolute ascending i.e smaller first).
	for (i = 0; i < (n - 1); i++) {
		k = i;
		w = e[(n * k) + k];
		for (j = (i + 1); j < n; j++) {
			if (mc_fabs(e[(n * j) + j]) < mc_fabs(w)) {
				k = j;
				w = e[(n * k) + k];
			}
		}
		if (k > i) {
			e[(n * k) + k] = e[(n * i) + i];
			e[(n * i) + i] = w;
			if (wantv) {
				for (j = 0; j < n; j++) {
					w              = v[(n * j) + k];
					v[(n * j) + k] = v[(n * j) + i];
					v[(n * j) + i] = w;
				}
			}
		}
	}

	return 0;
}

MC_TARGET_FUNC int mc_jacobisynxnl(const int n, const long double * a, long double tol, long double * e, long double * v)
{
//!# Requires a[n x n], e[n x n] and v[n x n] if !null where 1 < n.
//!# A and E may be the same. Computing eigenvalues and right-hand
//!# side eigenvectors of a real symmetric matrix by Jacobi method.
	int i, j, k, l = 0;
	long double c, s, scale, sumsq, diff, w, t;

	const int max   = 190;
	const int wantv = mc_nonnullptr(v);

	if (a != e) {
		mc_copynxnl(n, e, a);
	}

	if (wantv) {
		mc_eyenxnl(n, v, 0);
	}

	if (tol <= 0.0L) {
		tol = mc_fabsl(tol);
	}
	if (tol == 0.0L) {
		tol = MCLIMITS_TINYL;
	}

	while (1) {
//!# Computing the partial Frobenius norm, iterating the
//!# upper-triangle excluding the main-diagonal elements.
		mc_triussqrnxnl(n, e, &sumsq, &scale, 1);
		if (mc_fisnearl(scale * mc_sqrtl(sumsq * 2.0L), 0.0L, 1)) {
			break;
		}
		for (i = 0; i < n; i++) {
			for (j = i + 1; j < n; j++) {
				t = e[(n * i) + j];
				c = 1.0L;
				s = 0.0L;
//!# Forming symmetric Schur.
				if (mc_fabsl(t) > tol) {
					diff = (e[(n * j) + j] - e[(n * i) + i]) / (2.0L * t);
					t    = (mc_copysignl(1.0L, diff) > 0.0L
						? ( 1.0L / ( diff + mc_hypot2l(1.0L, diff)))
						: (-1.0L / (-diff + mc_hypot2l(1.0L, diff)))
					);
					w = 1.0L / mc_hypot2l(1.0L, t);
					c = w;
					s = t * w;
				}
//!# Rotating left.
				for (k = 0; k < n; k++) {
					t              = e[(n * i) + k];
					w              = e[(n * j) + k];
					e[(n * i) + k] = t * c - w * s;
					e[(n * j) + k] = t * s + w * c;
				}
//!# Rotating right.
				for (k = 0; k < n; k++) {
					t              = e[(n * k) + i];
					w              = e[(n * k) + j];
					e[(n * k) + i] = t * c - w * s;
					e[(n * k) + j] = t * s + w * c;
					if (wantv) {
						t              = v[(n * k) + i];
						w              = v[(n * k) + j];
						v[(n * k) + i] = t * c - w * s;
						v[(n * k) + j] = t * s + w * c;
					}
				}
			}
		}
		++l;
		if (l >= max) {
			return -1;
		}
	}

//!# Reordering eigenvalues and eigenvectors (absolute ascending i.e smaller first).
	for (i = 0; i < (n - 1); i++) {
		k = i;
		w = e[(n * k) + k];
		for (j = (i + 1); j < n; j++) {
			if (mc_fabsl(e[(n * j) + j]) < mc_fabsl(w)) {
				k = j;
				w = e[(n * k) + k];
			}
		}
		if (k > i) {
			e[(n * k) + k] = e[(n * i) + i];
			e[(n * i) + i] = w;
			if (wantv) {
				for (j = 0; j < n; j++) {
					w              = v[(n * j) + k];
					v[(n * j) + k] = v[(n * j) + i];
					v[(n * j) + i] = w;
				}
			}
		}
	}

	return 0;
}

#endif /* !MC_JACOBISYNXN_H */

/* EOF */
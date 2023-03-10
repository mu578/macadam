//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_svdgr1mxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_hypot2.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/mcswap.h>

#ifndef MC_SVDGRMXN_H
#define MC_SVDGRMXN_H

#pragma mark - mc_svdgr1mxn -

MC_TARGET_FUNC int mc_svdgr1mxnf(const int m, const int n, const float * a, float * w, int withu, int withv, int sorted, float eps, float tol, float * u, float * MC_TARGET_RESTRICT s, float * MC_TARGET_RESTRICT v)
{
//!# A and U may be the same. Requires a[m x n], w[n], u[m x p], s[1 x p] and v[p x p] where 1 < n <= m hence p=n.
//!# G. H. GOLUB and C. REINSCH SVD algorithm by bidiagonalization and modified QR steps. Handbook for Automatic
//!# Computation, vol. II, Linear Algebra", Springer-Verlag. Uses Householder transformations to reduce A to
//!# bidiagonal form, and then the QR algorithm to find the singular values of the bidiagonal matrix. The two
//!# phases properly combined produce the singular value decomposition of A.
//!#
//!# The main result SVD provides is that we can write an m by n matrix A
//!# such as U'*A=S*V' with:
//!#     - U is an [m x p] orthogonal matrix. The left-singular vectors of A are a set of orthonormal eigenvectors of AA'.
//!#     - S is an [n x p] diagonal matrix. The non-negative singular values of A (found on the diagonal entries of S) are
//!#       the square roots of the non-negative eigenvalues of both AA' and A'A.
//!#     - V is an [p x p] orthogonal matrix. The right-singular vectors of A are a set of orthonormal eigenvectors of A'A.
//!#     - p=min(m, n) and in this particular case we have m >= n then p=n; S is stored as a vector diagonal of size [1 x p].
//!#
//!# \note This algorithm would statify 95% of regular people use cases, however it
//!# won't compete, in term of error, against a double divide-and-conquer approach.
//!#
//!# Relative desired machine float-min. Optimal tuning: tol = TINY_VALUE greater than true FMIN.
//!# Relative desired machine epsilon.   Optimal tuning: eps = 2 * MACHINE_EPSILON.
//!#

//!# Max iteration allowed for convergence.
	const int max = 90;

	int i, j, k, l, q, iter, r = -1;
	float c, e, f, g = 0.0f, h, x, y, z;

	if (m > 1 && m >= n) {
		r = 0;
		if (a != u) {
			for (i = 0; i < (m * n); i++) {
				if (i < n) {
					w[i] = 0.0f;
				}
				u[i] = a[i];
			}
		}
//!# Step 1: Householder's reduction to bidiagonal form.
		x = 0.0f;
		for (i = 0; i < n; i++) {
			w[i] = g;
			e    = 0.0f;
			l    = i + 1;
			for (j = i; j < m; j++) {
				e = e + (u[(n * j) + i] * u[(n * j) + i]);
			}
			if (e < tol) {
				g = 0.0f;
			} else {
				f              =  u[(n * i) + i];
				g              = -mc_copysignf(1.0f, f) * mc_sqrtf(e);
				h              =  f * g - e;
				u[(n * i) + i] =  f - g;
				for (j = l; j < n; j++) {
					e = 0.0f;
					for (k = i; k < m; k++) {
						e = e + (u[(n * k) + i] * u[(n * k) + j]);
					}
					if (h == 0.0f) {
						continue;
					}
					f = e / h;
					for (k = i; k < m; k++) {
						u[(n * k) + j] = u[(n * k) + j] + (f * u[(n * k) + i]);
					}
				}
			}
			s[i] = g;
			e    = 0.0f;
			for (j = l; j < n; j++) {
				e = e + (u[(n * i) + j] * u[(n * i) + j]);
			}
			if (e < tol) {
				g = 0.0f;
			} else {
				f =  u[(n * i) + (i + 1)];
				g = -mc_copysignf(1.0f, f) * mc_sqrtf(e);
				h =  f * g - e;
				if (h == 0.0f) {
					continue;
				}
				u[(n * i) + (i + 1)] = f - g;
				for (j = l; j < n; j++) {
					w[j] = u[(n * i) + j] / h;
				}
				for (j = l; j < m; j++) {
					e = 0.0f;
					for (k = l; k < n; k++) {
						e = e + (u[(n * j) + k] * u[(n * i) + k]);
					}
					for (k = l; k < n; k++) {
						u[(n * j) + k] = u[(n * j) + k] + (e * w[k]);
					}
				}
			}
			y = mc_fabsf(s[i]) + mc_fabsf(w[i]);
			if (y > x) {
				x = y;
			}
		}

//!# Step 2: Accumulation of right-hand transformations.
		if (withv) {
			for (i = (n - 1); i >= 0; i--) {
				if (g != 0.0f) {
					h = u[(n * i) + (i + 1)] * g;
					if (h == 0.0f) {
						continue;
					}
					for (j = l; j < n; j++) {
						v[(n * j) + i] = u[(n * i) + j] / h;
					}
					for (j = l; j < n; j++) {
						e = 0.0f;
						for (k = l; k < n; k++) {
							e = e + (u[(n * i) + k] * v[(n * k) + j]);
						}
						for (k = l; k < n; k++) {
							v[(n * k) + j] = v[(n * k) + j] + (e * v[(n * k) + i]);
						}
					}
				}
				for (j = l; j < n; j++) {
					v[(n * i) + j] = v[(n * j) + i] = 0.0f;
				}
				v[(n * i) + i] = 1.0f;
				g              = w[i];
				l              = i;
			}
		}

//!# Step 3: Accumulation of left-hand transformations.
		if (withu) {
			for (i = n - 1; i >= 0; --i) {
				l = i + 1;
				g = s[i];
				if (i < (n - 1)) {
					for (j = l; j < n; ++j) {
						u[(n * i) + j] = 0.0f;
					}
				}
				if (g != 0.0f) {
					g = 1.0f / g;
					if (i != (n - 1)) {
						for (j = l; j < n; ++j) {
							h = 0.0f;
							for (k = l; k < m; ++k) {
								h = h + u[(n * k) + i] * u[(n * k) + j];
							}
							f = (h / u[(n * i) + i]) * g;
							for (k = i; k < m; ++k) {
								u[(n * k) + j] = u[(n * k) + j] + f * u[(n * k) + i];
							}
						}
					}
					for (j = i; j < m; ++j) {
						u[(n * j) + i] = u[(n * j) + i] * g;
					}
				} else {
					for (j = i; j < m; ++j) {
						u[(n * j) + i] = 0.0f;
					}
				}
				u[(n * i) + i] = u[(n * i) + i] + 1.0f;
			}
		}

//!# Step 4: Diagonalization of the bidiagonal form, extracting singular-values
//!#         and finalizing U such as U=A*V*S-1 by performing QR-steps.
		eps = eps * x;
		for (k = (n - 1); k >= 0; k--) {
			iter = 0;
	lnext:
//!# Testing convergence or cancellation.
			for (l = k; l >= 0; l--) {
				if (mc_fabsf(w[l]) <= eps) {
					goto l20;
				}
				if (mc_fabsf(s[l - 1]) <= eps) {
					goto l10;
				}
			}
	l10:
//!# Cancellation of e[l] if l > 0.
			c = 0.0f;
			e = 1.0f;
			q = l - 1;
			for (i = l; i <= k; i++) {
				f    = e * w[i];
				w[i] = w[i] * c;
				if (mc_fabsf(f) <= eps) {
					goto l20;
				}
				g  = s[i];
				h  = s[i] = mc_hypot2f(f, g);
				if (h == 0.0f) {
					continue;
				}
				c  = g / h;
				e = -f / h;
				if (withu) {
					for (j = 0; j < m; j++) {
						y              = u[(n * j) + q];
						z              = u[(n * j) + i];
						u[(n * j) + q] =  y * c + z * e;
						u[(n * j) + i] = -y * e + z * c;
					}
				}
			}
	l20:
//!# Testing convergence or exiting.
			z = s[k];
			if (l == k) {
				goto l30;
			}
//!# Enforcing hardlimit on convergence iterations.
			++iter;
			if (iter > max) {
				r = -1;
				goto l40;
				break;
			}

//!# Shift from bottom 2x2 minors.
			x = s[l];
			y = s[k - 1];
			g = w[k - 1];
			h = w[k];
			if (y == 0.0f || h == 0.0f) {
				continue;
			}
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0f * h * y);
			g = mc_hypot2f(f, 1.0f);
			f = ((x - z) * (x + z) + h * (y / ((f < 0) ? (f - g) : (f + g)) - h)) / x;

//!# Francis QR step.
			c = 1.0f;
			e = 1.0f;
			for (i = l + 1; i <= k;i++) {
				g        = w[i];
				y        = s[i];
				h        = e * g;
				g        = g * c;
				z        = mc_hypot2f(f, h);
				if (z == 0.0f) {
					continue;
				}
				w[i - 1] =  z;
				c        =  f / z;
				e        =  h / z;
				f        =  x * c + g * e;
				g        = -x * e + g * c;
				h        =  y * e;
				y        =  y * c;
				if (withv) {
					for (j = 0; j < n; j++) {
						x                    =  v[(n * j) + (i - 1)];
						z                    =  v[(n * j) + i];
						v[(n * j) + (i - 1)] =  x * c + z * e;
						v[(n * j) + i]       = -x * e + z * c;
					}
				}
				z        = mc_hypot2f(f, h);
				if (z == 0.0f) {
					continue;
				}
				s[i - 1] =  z;
				c        =  f / z;
				e        =  h / z;
				f        =  c * g + e * y;
				x        = -e * g + c * y;
				if (withu) {
					for (j = 0; j < m; j++) {
						y                    =  u[(n * j) + (i - 1)];
						z                    =  u[(n * j) + i];
						u[(n * j) + (i - 1)] =  y * c + z * e;
						u[(n * j) + i]       = -y * e + z * c;
					}
				}
			}
			w[l] = 0.0f;
			w[k] = f;
			s[k] = x;
//!# Performs next step if any.
			goto lnext;
	l30:
//!# Convergence, changing sign if necessary.
			if (mc_copysignf(1.0f, z) < 0.0f) {
				s[k] = -z;
				if (withv) {
					for (j = 0; j < n; j++) {
						v[(n * j) + k] = -v[(n * j) + k];
					}
				}
			}
		}
		if (sorted) {
//!# Step 6: Reordering descending.
			for (i = 0; i < n - 1; i++) {
				for (j = 0; j < n - 1 - i; j++) {
					if (s[j] < s[j + 1]) {
						mcswap_var(x, s[j], s[j + 1]);
						for (k = 0; k < m; k++) {
							if (withv && k < n) {
								mcswap_var(x, v[(n * k) + j], v[(n * k) + j + 1]);
							}
							if (withu) {
								mcswap_var(x, u[(n * k) + j], u[(n * k) + j + 1]);
							} else if (k >= (n - 1)) {
								break;
							}
						}
					}
				}
			}
		}
	}
l40:
	return r;
}

MC_TARGET_FUNC int mc_svdgr1mxnff(const int m, const int n, const float * a, double * w, int withu, int withv, int sorted, double eps, double tol, double * u, double * MC_TARGET_RESTRICT s, double * MC_TARGET_RESTRICT v)
{
//!# A and U may be the same. Requires a[m x n], w[n], u[m x p], s[1 x p] and v[p x p] where 1 < n <= m hence p=n.
//!# G. H. GOLUB and C. REINSCH SVD algorithm by bidiagonalization and modified QR steps. Handbook for Automatic
//!# Computation, vol. II, Linear Algebra", Springer-Verlag. Uses Householder transformations to reduce A to
//!# bidiagonal form, and then the QR algorithm to find the singular values of the bidiagonal matrix. The two
//!# phases properly combined produce the singular value decomposition of A.
//!#
//!# The main result SVD provides is that we can write an m by n matrix A
//!# such as U'*A=S*V' with:
//!#     - U is an [m x p] orthogonal matrix. The left-singular vectors of A are a set of orthonormal eigenvectors of AA'.
//!#     - S is an [n x p] diagonal matrix. The non-negative singular values of A (found on the diagonal entries of S) are
//!#       the square roots of the non-negative eigenvalues of both AA' and A'A.
//!#     - V is an [p x p] orthogonal matrix. The right-singular vectors of A are a set of orthonormal eigenvectors of A'A.
//!#     - p=min(m, n) and in this particular case we have m >= n then p=n; S is stored as a vector diagonal of size [1 x p].
//!#
//!# \note This algorithm would statify 95% of regular people use cases, however it
//!# won't compete, in term of error, against a double divide-and-conquer approach.
//!#
//!# Relative desired machine double-min. Optimal tuning: tol = TINY_VALUE greater than true FMIN.
//!# Relative desired machine epsilon.   Optimal tuning: eps = 2 * MACHINE_EPSILON.
//!#

//!# Max iteration allowed for convergence.
	const int max = 120;

	int i, j, k, l, q, iter, r = -1;
	double c, e, f, g = 0.0, h, x, y, z;

	if (m > 1 && m >= n) {
		r = 0;
		for (i = 0; i < (m * n); i++) {
			if (i < n) {
				w[i] = 0.0;
			}
			u[i] = mc_cast(double, a[i]);
		}
//!# Step 1: Householder's reduction to bidiagonal form.
		x = 0.0;
		for (i = 0; i < n; i++) {
			w[i] = g;
			e    = 0.0;
			l    = i + 1;
			for (j = i; j < m; j++) {
				e = e + (u[(n * j) + i] * u[(n * j) + i]);
			}
			if (e < tol) {
				g = 0.0;
			} else {
				f              =  u[(n * i) + i];
				g              = -mc_copysign(1.0, f) * mc_sqrt(e);
				h              =  f * g - e;
				u[(n * i) + i] =  f - g;
				for (j = l; j < n; j++) {
					e = 0.0;
					for (k = i; k < m; k++) {
						e = e + (u[(n * k) + i] * u[(n * k) + j]);
					}
					if (h == 0.0) {
						continue;
					}
					f = e / h;
					for (k = i; k < m; k++) {
						u[(n * k) + j] = u[(n * k) + j] + (f * u[(n * k) + i]);
					}
				}
			}
			s[i] = g;
			e    = 0.0;
			for (j = l; j < n; j++) {
				e = e + (u[(n * i) + j] * u[(n * i) + j]);
			}
			if (e < tol) {
				g = 0.0;
			} else {
				f =  u[(n * i) + (i + 1)];
				g = -mc_copysign(1.0, f) * mc_sqrt(e);
				h =  f * g - e;
				if (h == 0.0) {
					continue;
				}
				u[(n * i) + (i + 1)] = f - g;
				for (j = l; j < n; j++) {
					w[j] = u[(n * i) + j] / h;
				}
				for (j = l; j < m; j++) {
					e = 0.0;
					for (k = l; k < n; k++) {
						e = e + (u[(n * j) + k] * u[(n * i) + k]);
					}
					for (k = l; k < n; k++) {
						u[(n * j) + k] = u[(n * j) + k] + (e * w[k]);
					}
				}
			}
			y = mc_fabs(s[i]) + mc_fabs(w[i]);
			if (y > x) {
				x = y;
			}
		}

//!# Step 2: Accumulation of right-hand transformations.
		if (withv) {
			for (i = (n - 1); i >= 0; i--) {
				if (g != 0.0) {
					h = u[(n * i) + (i + 1)] * g;
					if (h == 0.0) {
						continue;
					}
					for (j = l; j < n; j++) {
						v[(n * j) + i] = u[(n * i) + j] / h;
					}
					for (j = l; j < n; j++) {
						e = 0.0;
						for (k = l; k < n; k++) {
							e = e + (u[(n * i) + k] * v[(n * k) + j]);
						}
						for (k = l; k < n; k++) {
							v[(n * k) + j] = v[(n * k) + j] + (e * v[(n * k) + i]);
						}
					}
				}
				for (j = l; j < n; j++) {
					v[(n * i) + j] = v[(n * j) + i] = 0.0;
				}
				v[(n * i) + i] = 1.0;
				g              = w[i];
				l              = i;
			}
		}

//!# Step 3: Accumulation of left-hand transformations.
		if (withu) {
			for (i = n - 1; i >= 0; --i) {
				l = i + 1;
				g = s[i];
				if (i < (n - 1)) {
					for (j = l; j < n; ++j) {
						u[(n * i) + j] = 0.0;
					}
				}
				if (g != 0.0) {
					g = 1.0 / g;
					if (i != (n - 1)) {
						for (j = l; j < n; ++j) {
							h = 0.0;
							for (k = l; k < m; ++k) {
								h = h + u[(n * k) + i] * u[(n * k) + j];
							}
							f = (h / u[(n * i) + i]) * g;
							for (k = i; k < m; ++k) {
								u[(n * k) + j] = u[(n * k) + j] + f * u[(n * k) + i];
							}
						}
					}
					for (j = i; j < m; ++j) {
						u[(n * j) + i] = u[(n * j) + i] * g;
					}
				} else {
					for (j = i; j < m; ++j) {
						u[(n * j) + i] = 0.0;
					}
				}
				u[(n * i) + i] = u[(n * i) + i] + 1.0;
			}
		}

//!# Step 4: Diagonalization of the bidiagonal form, extracting singular-values
//!#         and finalizing U such as U=A*V*S-1 by performing QR-steps.
		eps = eps * x;
		for (k = (n - 1); k >= 0; k--) {
			iter = 0;
	lnext:
//!# Testing convergence or cancellation.
			for (l = k; l >= 0; l--) {
				if (mc_fabs(w[l]) <= eps) {
					goto l20;
				}
				if (mc_fabs(s[l - 1]) <= eps) {
					goto l10;
				}
			}
	l10:
//!# Cancellation of e[l] if l > 0.
			c = 0.0;
			e = 1.0;
			q = l - 1;
			for (i = l; i <= k; i++) {
				f    = e * w[i];
				w[i] = w[i] * c;
				if (mc_fabs(f) <= eps) {
					goto l20;
				}
				g  = s[i];
				h  = s[i] = mc_hypot2(f, g);
				if (h == 0.0) {
					continue;
				}
				c  = g / h;
				e = -f / h;
				if (withu) {
					for (j = 0; j < m; j++) {
						y              =  u[(n * j) + q];
						z              =  u[(n * j) + i];
						u[(n * j) + q] =  y * c + z * e;
						u[(n * j) + i] = -y * e + z * c;
					}
				}
			}
	l20:
//!# Testing convergence or exiting.
			z = s[k];
			if (l == k) {
				goto l30;
			}
//!# Enforcing hardlimit on convergence iterations.
			++iter;
			if (iter > max) {
				r = -1;
				goto l40;
				break;
			}

//!# Shift from bottom 2x2 minors.
			x = s[l];
			y = s[k - 1];
			g = w[k - 1];
			h = w[k];
			if (y == 0.0 || h == 0.0) {
				continue;
			}
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = mc_hypot2(f, 1.0);
			f = ((x - z) * (x + z) + h * (y / ((f < 0) ? (f - g) : (f + g)) - h)) / x;

//!# Francis QR step.
			c = 1.0;
			e = 1.0;
			for (i = l + 1; i <= k;i++) {
				g        = w[i];
				y        = s[i];
				h        = e * g;
				g        = g * c;
				z        = mc_hypot2(f, h);
				if (z == 0.0) {
					continue;
				}
				w[i - 1] =  z;
				c        =  f / z;
				e        =  h / z;
				f        =  x * c + g * e;
				g        = -x * e + g * c;
				h        =  y * e;
				y        =  y * c;
				if (withv) {
					for (j = 0; j < n; j++) {
						x                    =  v[(n * j) + (i - 1)];
						z                    =  v[(n * j) + i];
						v[(n * j) + (i - 1)] =  x * c + z * e;
						v[(n * j) + i]       = -x * e + z * c;
					}
				}
				z        = mc_hypot2(f, h);
				if (z == 0.0) {
					continue;
				}
				s[i - 1] =  z;
				c        =  f / z;
				e        =  h / z;
				f        =  c * g + e * y;
				x        = -e * g + c * y;
				if (withu) {
					for (j = 0; j < m; j++) {
						y                    =  u[(n * j) + (i - 1)];
						z                    =  u[(n * j) + i];
						u[(n * j) + (i - 1)] =  y * c + z * e;
						u[(n * j) + i]       = -y * e + z * c;
					}
				}
			}
			w[l] = 0.0;
			w[k] = f;
			s[k] = x;
//!# Performs next step if any.
			goto lnext;
	l30:
//!# Convergence, changing sign if necessary.
			if (mc_copysign(1.0, z) < 0.0) {
				s[k] = -z;
				if (withv) {
					for (j = 0; j < n; j++) {
						v[(n * j) + k] = -v[(n * j) + k];
					}
				}
			}
		}
		if (sorted) {
//!# Step 6: Reordering descending.
			for (i = 0; i < n - 1; i++) {
				for (j = 0; j < n - 1 - i; j++) {
					if (s[j] < s[j + 1]) {
						mcswap_var(x, s[j], s[j + 1]);
						for (k = 0; k < m; k++) {
							if (withv && k < n) {
								mcswap_var(x, v[(n * k) + j], v[(n * k) + j + 1]);
							}
							if (withu) {
								mcswap_var(x, u[(n * k) + j], u[(n * k) + j + 1]);
							} else if (k >= (n - 1)) {
								break;
							}
						}
					}
				}
			}
		}
	}
l40:
	return r;
}

MC_TARGET_FUNC int mc_svdgr1mxn(const int m, const int n, const double * a, double * w, int withu, int withv, int sorted, double eps, double tol, double * u, double * MC_TARGET_RESTRICT s, double * MC_TARGET_RESTRICT v)
{
//!# A and U may be the same. Requires a[m x n], w[n], u[m x p], s[1 x p] and v[p x p] where 1 < n <= m hence p=n.
//!# G. H. GOLUB and C. REINSCH SVD algorithm by bidiagonalization and modified QR steps. Handbook for Automatic
//!# Computation, vol. II, Linear Algebra", Springer-Verlag. Uses Householder transformations to reduce A to
//!# bidiagonal form, and then the QR algorithm to find the singular values of the bidiagonal matrix. The two
//!# phases properly combined produce the singular value decomposition of A.
//!#
//!# The main result SVD provides is that we can write an m by n matrix A
//!# such as U'*A=S*V' with:
//!#     - U is an [m x p] orthogonal matrix. The left-singular vectors of A are a set of orthonormal eigenvectors of AA'.
//!#     - S is an [n x p] diagonal matrix. The non-negative singular values of A (found on the diagonal entries of S) are
//!#       the square roots of the non-negative eigenvalues of both AA' and A'A.
//!#     - V is an [p x p] orthogonal matrix. The right-singular vectors of A are a set of orthonormal eigenvectors of A'A.
//!#     - p=min(m, n) and in this particular case we have m >= n then p=n; S is stored as a vector diagonal of size [1 x p].
//!#
//!# \note This algorithm would statify 95% of regular people use cases, however it
//!# won't compete, in term of error, against a double divide-and-conquer approach.
//!#
//!# Relative desired machine double-min. Optimal tuning: tol = TINY_VALUE greater than true FMIN.
//!# Relative desired machine epsilon.   Optimal tuning: eps = 2 * MACHINE_EPSILON.
//!#

//!# Max iteration allowed for convergence.
	const int max = 120;

	int i, j, k, l, q, iter, r = -1;
	double c, e, f, g = 0.0, h, x, y, z;

	if (m > 1 && m >= n) {
		r = 0;
		if (a != u) {
			for (i = 0; i < (m * n); i++) {
				if (i < n) {
					w[i] = 0.0;
				}
				u[i] = a[i];
			}
		}
//!# Step 1: Householder's reduction to bidiagonal form.
		x = 0.0;
		for (i = 0; i < n; i++) {
			w[i] = g;
			e    = 0.0;
			l    = i + 1;
			for (j = i; j < m; j++) {
				e = e + (u[(n * j) + i] * u[(n * j) + i]);
			}
			if (e < tol) {
				g = 0.0;
			} else {
				f              =  u[(n * i) + i];
				g              = -mc_copysign(1.0, f) * mc_sqrt(e);
				h              =  f * g - e;
				u[(n * i) + i] =  f - g;
				for (j = l; j < n; j++) {
					e = 0.0;
					for (k = i; k < m; k++) {
						e = e + (u[(n * k) + i] * u[(n * k) + j]);
					}
					if (h == 0.0) {
						continue;
					}
					f = e / h;
					for (k = i; k < m; k++) {
						u[(n * k) + j] = u[(n * k) + j] + (f * u[(n * k) + i]);
					}
				}
			}
			s[i] = g;
			e    = 0.0;
			for (j = l; j < n; j++) {
				e = e + (u[(n * i) + j] * u[(n * i) + j]);
			}
			if (e < tol) {
				g = 0.0;
			} else {
				f =  u[(n * i) + (i + 1)];
				g = -mc_copysign(1.0, f) * mc_sqrt(e);
				h =  f * g - e;
				if (h == 0.0) {
					continue;
				}
				u[(n * i) + (i + 1)] = f - g;
				for (j = l; j < n; j++) {
					w[j] = u[(n * i) + j] / h;
				}
				for (j = l; j < m; j++) {
					e = 0.0;
					for (k = l; k < n; k++) {
						e = e + (u[(n * j) + k] * u[(n * i) + k]);
					}
					for (k = l; k < n; k++) {
						u[(n * j) + k] = u[(n * j) + k] + (e * w[k]);
					}
				}
			}
			y = mc_fabs(s[i]) + mc_fabs(w[i]);
			if (y > x) {
				x = y;
			}
		}

//!# Step 2: Accumulation of right-hand transformations.
		if (withv) {
			for (i = (n - 1); i >= 0; i--) {
				if (g != 0.0) {
					h = u[(n * i) + (i + 1)] * g;
					if (h == 0.0) {
						continue;
					}
					for (j = l; j < n; j++) {
						v[(n * j) + i] = u[(n * i) + j] / h;
					}
					for (j = l; j < n; j++) {
						e = 0.0;
						for (k = l; k < n; k++) {
							e = e + (u[(n * i) + k] * v[(n * k) + j]);
						}
						for (k = l; k < n; k++) {
							v[(n * k) + j] = v[(n * k) + j] + (e * v[(n * k) + i]);
						}
					}
				}
				for (j = l; j < n; j++) {
					v[(n * i) + j] = v[(n * j) + i] = 0.0;
				}
				v[(n * i) + i] = 1.0;
				g              = w[i];
				l              = i;
			}
		}

//!# Step 3: Accumulation of left-hand transformations.
		if (withu) {
			for (i = n - 1; i >= 0; --i) {
				l = i + 1;
				g = s[i];
				if (i < (n - 1)) {
					for (j = l; j < n; ++j) {
						u[(n * i) + j] = 0.0;
					}
				}
				if (g != 0.0) {
					g = 1.0 / g;
					if (i != (n - 1)) {
						for (j = l; j < n; ++j) {
							h = 0.0;
							for (k = l; k < m; ++k) {
								h = h + u[(n * k) + i] * u[(n * k) + j];
							}
							f = (h / u[(n * i) + i]) * g;
							for (k = i; k < m; ++k) {
								u[(n * k) + j] = u[(n * k) + j] + f * u[(n * k) + i];
							}
						}
					}
					for (j = i; j < m; ++j) {
						u[(n * j) + i] = u[(n * j) + i] * g;
					}
				} else {
					for (j = i; j < m; ++j) {
						u[(n * j) + i] = 0.0;
					}
				}
				u[(n * i) + i] = u[(n * i) + i] + 1.0;
			}
		}

//!# Step 4: Diagonalization of the bidiagonal form, extracting singular-values
//!#         and finalizing U such as U=A*V*S-1 by performing QR-steps.
		eps = eps * x;
		for (k = (n - 1); k >= 0; k--) {
			iter = 0;
	lnext:
//!# Testing convergence or cancellation.
			for (l = k; l >= 0; l--) {
				if (mc_fabs(w[l]) <= eps) {
					goto l20;
				}
				if (mc_fabs(s[l - 1]) <= eps) {
					goto l10;
				}
			}
	l10:
//!# Cancellation of e[l] if l > 0.
			c = 0.0;
			e = 1.0;
			q = l - 1;
			for (i = l; i <= k; i++) {
				f    = e * w[i];
				w[i] = w[i] * c;
				if (mc_fabs(f) <= eps) {
					goto l20;
				}
				g  = s[i];
				h  = s[i] = mc_hypot2(f, g);
				if (h == 0.0) {
					continue;
				}
				c  = g / h;
				e = -f / h;
				if (withu) {
					for (j = 0; j < m; j++) {
						y              =  u[(n * j) + q];
						z              =  u[(n * j) + i];
						u[(n * j) + q] =  y * c + z * e;
						u[(n * j) + i] = -y * e + z * c;
					}
				}
			}
	l20:
//!# Testing convergence or exiting.
			z = s[k];
			if (l == k) {
				goto l30;
			}
//!# Enforcing hardlimit on convergence iterations.
			++iter;
			if (iter > max) {
				r = -1;
				goto l40;
				break;
			}

//!# Shift from bottom 2x2 minors.
			x = s[l];
			y = s[k - 1];
			g = w[k - 1];
			h = w[k];
			if (mc_fabs(y) == 0.0 || h == 0.0) {
				continue;
			}
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = mc_hypot2(f, 1.0);
			f = ((x - z) * (x + z) + h * (y / ((f < 0) ? (f - g) : (f + g)) - h)) / x;

//!# Francis QR step.
			c = 1.0;
			e = 1.0;
			for (i = l + 1; i <= k;i++) {
				g        = w[i];
				y        = s[i];
				h        = e * g;
				g        = g * c;
				z        = mc_hypot2(f, h);
				if (z == 0.0) {
					continue;
				}
				w[i - 1] =  z;
				c        =  f / z;
				e        =  h / z;
				f        =  x * c + g * e;
				g        = -x * e + g * c;
				h        =  y * e;
				y        =  y * c;
				if (withv) {
					for (j = 0; j < n; j++) {
						x                    =  v[(n * j) + (i - 1)];
						z                    =  v[(n * j) + i];
						v[(n * j) + (i - 1)] =  x * c + z * e;
						v[(n * j) + i]       = -x * e + z * c;
					}
				}
				z        = mc_hypot2(f, h);
				if (z == 0.0) {
					continue;
				}
				s[i - 1] = z;
				c        = f / z;
				e        = h / z;
				f        =  c * g + e * y;
				x        = -e * g + c * y;
				if (withu) {
					for (j = 0; j < m; j++) {
						y                    =  u[(n * j) + (i - 1)];
						z                    =  u[(n * j) + i];
						u[(n * j) + (i - 1)] =  y * c + z * e;
						u[(n * j) + i]       = -y * e + z * c;
					}
				}
			}
			w[l] = 0.0;
			w[k] = f;
			s[k] = x;
//!# Performs next step if any.
			goto lnext;
	l30:
//!# Convergence, changing sign if necessary.
			if (mc_copysign(1.0, z) < 0.0) {
				s[k] = -z;
				if (withv) {
					for (j = 0; j < n; j++) {
						v[(n * j) + k] = -v[(n * j) + k];
					}
				}
			}
		}
		if (sorted) {
//!# Step 6: Reordering descending.
			for (i = 0; i < n - 1; i++) {
				for (j = 0; j < n - 1 - i; j++) {
					if (s[j] < s[j + 1]) {
						mcswap_var(x, s[j], s[j + 1]);
						for (k = 0; k < m; k++) {
							if (withv && k < n) {
								mcswap_var(x, v[(n * k) + j], v[(n * k) + j + 1]);
							}
							if (withu) {
								mcswap_var(x, u[(n * k) + j], u[(n * k) + j + 1]);
							} else if (k >= (n - 1)) {
								break;
							}
						}
					}
				}
			}
		}
	}
l40:
	return r;
}

MC_TARGET_FUNC int mc_svdgr1mxnl(const int m, const int n, const long double * a, long double * w, int withu, int withv, int sorted, long double eps, long double tol, long double * u, long double * MC_TARGET_RESTRICT s, long double * MC_TARGET_RESTRICT v)
{
//!# A and U may be the same. Requires a[m x n], w[n], u[m x p], s[1 x p] and v[p x p] where 1 < n <= m hence p=n.
//!# G. H. GOLUB and C. REINSCH SVD algorithm by bidiagonalization and modified QR steps. Handbook for Automatic
//!# Computation, vol. II, Linear Algebra", Springer-Verlag. Uses Householder transformations to reduce A to
//!# bidiagonal form, and then the QR algorithm to find the singular values of the bidiagonal matrix. The two
//!# phases properly combined produce the singular value decomposition of A.
//!#
//!# The main result SVD provides is that we can write an m by n matrix A
//!# such as U'*A=S*V' with:
//!#     - U is an [m x p] orthogonal matrix. The left-singular vectors of A are a set of orthonormal eigenvectors of AA'.
//!#     - S is an [n x p] diagonal matrix. The non-negative singular values of A (found on the diagonal entries of S) are
//!#       the square roots of the non-negative eigenvalues of both AA' and A'A.
//!#     - V is an [p x p] orthogonal matrix. The right-singular vectors of A are a set of orthonormal eigenvectors of A'A.
//!#     - p=min(m, n) and in this particular case we have m >= n then p=n; S is stored as a vector diagonal of size [1 x p].
//!#
//!# \note This algorithm would statify 95% of regular people use cases, however it
//!# won't compete, in term of error, against a double divide-and-conquer approach.
//!#
//!# Relative desired machine long double-min. Optimal tuning: tol = TINY_VALUE greater than true FMIN.
//!# Relative desired machine epsilon.   Optimal tuning: eps = 2 * MACHINE_EPSILON.
//!#

//!# Max iteration allowed for convergence.
	const int max = 180;

	int i, j, k, l, q, iter, r = -1;
	long double c, e, f, g = 0.0L, h, x, y, z;

	if (m > 1 && m >= n) {
		r = 0;
		if (a != u) {
			for (i = 0; i < (m * n); i++) {
				if (i < n) {
					w[i] = 0.0L;
				}
				u[i] = a[i];
			}
		}
//!# Step 1: Householder's reduction to bidiagonal form.
		x = 0.0L;
		for (i = 0; i < n; i++) {
			w[i] = g;
			e    = 0.0L;
			l    = i + 1;
			for (j = i; j < m; j++) {
				e = e + (u[(n * j) + i] * u[(n * j) + i]);
			}
			if (e < tol) {
				g = 0.0L;
			} else {
				f              =  u[(n * i) + i];
				g              = -mc_copysignl(1.0L, f) * mc_sqrtl(e);
				h              =  f * g - e;
				u[(n * i) + i] =  f - g;
				for (j = l; j < n; j++) {
					e = 0.0L;
					for (k = i; k < m; k++) {
						e = e + (u[(n * k) + i] * u[(n * k) + j]);
					}
					if (h == 0.0L) {
						continue;
					}
					f = e / h;
					for (k = i; k < m; k++) {
						u[(n * k) + j] = u[(n * k) + j] + (f * u[(n * k) + i]);
					}
				}
			}
			s[i] = g;
			e    = 0.0L;
			for (j = l; j < n; j++) {
				e = e + (u[(n * i) + j] * u[(n * i) + j]);
			}
			if (e < tol) {
				g = 0.0L;
			} else {
				f =  u[(n * i) + (i + 1)];
				g = -mc_copysignl(1.0L, f) * mc_sqrtl(e);
				h =  f * g - e;
				if (h == 0.0L) {
					continue;
				}
				u[(n * i) + (i + 1)] = f - g;
				for (j = l; j < n; j++) {
					w[j] = u[(n * i) + j] / h;
				}
				for (j = l; j < m; j++) {
					e = 0.0L;
					for (k = l; k < n; k++) {
						e = e + (u[(n * j) + k] * u[(n * i) + k]);
					}
					for (k = l; k < n; k++) {
						u[(n * j) + k] = u[(n * j) + k] + (e * w[k]);
					}
				}
			}
			y = mc_fabsl(s[i]) + mc_fabsl(w[i]);
			if (y > x) {
				x = y;
			}
		}

//!# Step 2: Accumulation of right-hand transformations.
		if (withv) {
			for (i = (n - 1); i >= 0; i--) {
				if (g != 0.0L) {
					h = u[(n * i) + (i + 1)] * g;
					if (h == 0.0L) {
						continue;
					}
					for (j = l; j < n; j++) {
						v[(n * j) + i] = u[(n * i) + j] / h;
					}
					for (j = l; j < n; j++) {
						e = 0.0L;
						for (k = l; k < n; k++) {
							e = e + (u[(n * i) + k] * v[(n * k) + j]);
						}
						for (k = l; k < n; k++) {
							v[(n * k) + j] = v[(n * k) + j] + (e * v[(n * k) + i]);
						}
					}
				}
				for (j = l; j < n; j++) {
					v[(n * i) + j] = v[(n * j) + i] = 0.0L;
				}
				v[(n * i) + i] = 1.0L;
				g              = w[i];
				l              = i;
			}
		}

//!# Step 3: Accumulation of left-hand transformations.
		if (withu) {
			for (i = n - 1; i >= 0; --i) {
				l = i + 1;
				g = s[i];
				if (i < (n - 1)) {
					for (j = l; j < n; ++j) {
						u[(n * i) + j] = 0.0L;
					}
				}
				if (g != 0.0L) {
					g = 1.0L / g;
					if (i != (n - 1)) {
						for (j = l; j < n; ++j) {
							h = 0.0L;
							for (k = l; k < m; ++k) {
								h = h + u[(n * k) + i] * u[(n * k) + j];
							}
							f = (h / u[(n * i) + i]) * g;
							for (k = i; k < m; ++k) {
								u[(n * k) + j] = u[(n * k) + j] + f * u[(n * k) + i];
							}
						}
					}
					for (j = i; j < m; ++j) {
						u[(n * j) + i] = u[(n * j) + i] * g;
					}
				} else {
					for (j = i; j < m; ++j) {
						u[(n * j) + i] = 0.0L;
					}
				}
				u[(n * i) + i] = u[(n * i) + i] + 1.0L;
			}
		}

//!# Step 4: Diagonalization of the bidiagonal form, extracting singular-values
//!#         and finalizing U such as U=A*V*S-1 by performing QR-steps.
		eps = eps * x;
		for (k = (n - 1); k >= 0; k--) {
			iter = 0;
	lnext:
//!# Testing convergence or cancellation.
			for (l = k; l >= 0; l--) {
				if (mc_fabsl(w[l]) <= eps) {
					goto l20;
				}
				if (mc_fabsl(s[l - 1]) <= eps) {
					goto l10;
				}
			}
	l10:
//!# Cancellation of e[l] if l > 0.
			c = 0.0L;
			e = 1.0L;
			q = l - 1;
			for (i = l; i <= k; i++) {
				f    = e * w[i];
				w[i] = w[i] * c;
				if (mc_fabsl(f) <= eps) {
					goto l20;
				}
				g  = s[i];
				h  = s[i] = mc_hypot2l(f, g);
				if (h == 0.0L) {
					continue;
				}
				c  = g / h;
				e = -f / h;
				if (withu) {
					for (j = 0; j < m; j++) {
						y              =  u[(n * j) + q];
						z              =  u[(n * j) + i];
						u[(n * j) + q] =  y * c + z * e;
						u[(n * j) + i] = -y * e + z * c;
					}
				}
			}
	l20:
//!# Testing convergence or exiting.
			z = s[k];
			if (l == k) {
				goto l30;
			}
//!# Enforcing hardlimit on convergence iterations.
			++iter;
			if (iter > max) {
				r = -1;
				goto l40;
				break;
			}

//!# Shift from bottom 2x2 minors.
			x = s[l];
			y = s[k - 1];
			g = w[k - 1];
			h = w[k];
			if (y == 0.0L || h == 0.0L) {
				continue;
			}
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0L * h * y);
			g = mc_hypot2l(f, 1.0L);
			f = ((x - z) * (x + z) + h * (y / ((f < 0) ? (f - g) : (f + g)) - h)) / x;

//!# Francis QR step.
			c = 1.0L;
			e = 1.0L;
			for (i = l + 1; i <= k;i++) {
				g        = w[i];
				y        = s[i];
				h        = e * g;
				g        = g * c;
				z        = mc_hypot2l(f, h);
				if (z == 0.0L) {
					continue;
				}
				w[i - 1] =  z;
				c        =  f / z;
				e        =  h / z;
				f        =  x * c + g * e;
				g        = -x * e + g * c;
				h        =  y * e;
				y        =  y * c;
				if (withv) {
					for (j = 0; j < n; j++) {
						x                    =  v[(n * j) + (i - 1)];
						z                    =  v[(n * j) + i];
						v[(n * j) + (i - 1)] =  x * c + z * e;
						v[(n * j) + i]       = -x * e + z * c;
					}
				}
				z        = mc_hypot2l(f, h);
				if (z == 0.0L) {
					continue;
				}
				s[i - 1] = z;
				c        = f / z;
				e        = h / z;
				f        =  c * g + e * y;
				x        = -e * g + c * y;
				if (withu) {
					for (j = 0; j < m; j++) {
						y                    =  u[(n * j) + (i - 1)];
						z                    =  u[(n * j) + i];
						u[(n * j) + (i - 1)] =  y * c + z * e;
						u[(n * j) + i]       = -y * e + z * c;
					}
				}
			}
			w[l] = 0.0L;
			w[k] = f;
			s[k] = x;
//!# Performs next step if any.
			goto lnext;
	l30:
//!# Convergence, changing sign if necessary.
			if (mc_copysignl(1.0L, z) < 0.0L) {
				s[k] = -z;
				if (withv) {
					for (j = 0; j < n; j++) {
						v[(n * j) + k] = -v[(n * j) + k];
					}
				}
			}
		}
		if (sorted) {
//!# Step 6: Reordering descending.
			for (i = 0; i < n - 1; i++) {
				for (j = 0; j < n - 1 - i; j++) {
					if (s[j] < s[j + 1]) {
						mcswap_var(x, s[j], s[j + 1]);
						for (k = 0; k < m; k++) {
							if (withv && k < n) {
								mcswap_var(x, v[(n * k) + j], v[(n * k) + j + 1]);
							}
							if (withu) {
								mcswap_var(x, u[(n * k) + j], u[(n * k) + j + 1]);
							} else if (k >= (n - 1)) {
								break;
							}
						}
					}
				}
			}
		}
	}
l40:
	return r;
}

#endif /* !MC_SVDGRMXN_H */

/* EOF */
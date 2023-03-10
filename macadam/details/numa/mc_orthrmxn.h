//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_orthrmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fisnear.h>
#include <macadam/details/math/mc_fmax.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/numa/mc_copymxn.h>
#include <macadam/details/numa/mc_dotpmx1.h>
#include <macadam/details/numa/mc_eyenxn.h>
#include <macadam/details/numa/mc_l2normmx1.h>
#include <macadam/details/numa/mc_minmax1xn.h>
#include <macadam/details/numa/mc_zerosmx1.h>
#include <macadam/mcswap.h>

#ifndef MC_ORTHRMXN_H
#define MC_ORTHRMXN_H

#pragma mark - mc_orthrmxn -

MC_TARGET_FUNC int mc_orthrmxnf(const int m, const int n, const float * a, float tol, float * q, float * MC_TARGET_RESTRICT r, float * MC_TARGET_RESTRICT w, int * pvi)
{
//!# Requires a[m x n], q[m x n] and r[n x n] if !null and w[n] && pvi[n] if !null where 1 < n <= m.
//!# A and Q may be the same. Forming a ortho-normalized basis Q using Modified Gram-Schmidt method
//!# + a decimeting column step if norm < tol + iterative re-orthogonalization step for rank deficient
//!# systems. If R is not null upper-right-triangle is formed. @see Achiya Dax, `A modified Gram-schmidt
//!# algorithm with iterative orthogonalization and column pivoting`.
	const int wantr   = mc_nonnullptr(r);
	const int wantpvi = mc_nonnullptr(w) && mc_nonnullptr(pvi);

	int i, j, k, l, p;
	float bnorm, cnorm, dot, s;

	if (m >= n) {
		if (a != q) {
			mc_copymxnf(m, n, q, a);
		}
		if (wantr) {
			mc_eyenxnf(n, r, 0);
		}

		if (wantpvi) {
			for (j = 0; j < n; j++) {
				s      = mc_l2normmx1f(m, n, j, q);
				w[j]   = mc_raise2f(s);
				pvi[j] = j;
			}
		}

		if (tol <= 0.0f) {
			tol = mc_fabsf(tol);
		}
		if (tol == 0.0f) {
			tol = MCLIMITS_EPSILONF;
		}
		bnorm = 0.0f;
		for (k = 0; k < n; k++) {
//!# Step 1: pivoting if required.
			if (wantpvi) {
				mc_minmax1xnf(n - k, w + k, MC_NULLPTR, &s, MC_NULLPTR, &l);
				l = l + k;
				if (k != l) {
					for (i = 0; i < m; i++) {
						mcswap_var(s, q[(n * i) + k], q[(n * i) + l]);
						if (i < k) {
							mcswap_var(s, r[(n * i) + k], r[(n * i) + l]);
						}
					}
					mcswap_var(p, pvi[k], pvi[l]);
				}
			}
//!# Step 2: re-orthogonalization.
			if (k > 0) {
				for (i = 0; i < k; i++) {
					dot = mc_dotpmx1f(m, n, n, i, k, q, q, 1);
					if (wantr) {
						r[(n * i) + k] = r[(n * i) + k] + dot;
					}
					for (j = 0; j < m; j++) {
						q[(n * j) + k] = q[(n * j) + k] - (dot * q[(n * j) + i]);
					}
				}
			}
//!# Step 3: normalization.
			cnorm = mc_l2normmx1f(m, n, k, q);
			if (cnorm != 0.0f) {
//!# Step 4: close to zero decimation step.
				if (cnorm < tol * bnorm) {
					mc_zerosmx1f(m, n, k, q);
					q[(n * k) + k] = 1.0f;
					if (wantr) {
						r[(n * k) + k] = 0.0f;
					}
				} else {
					bnorm = mc_fmaxf(bnorm, cnorm);
					if (wantr) {
						r[(n * k) + k] = cnorm;
					}
					if (!mc_fisnearf(cnorm, 1.0f, 1)) {
						cnorm = 1.0f / cnorm;
						for (j = 0; j < m; j++) {
							q[(n * j) + k] = q[(n * j) + k] * cnorm;
						}
					}
				}
			} else {
				q[(n * k) + k] = 1.0f;
				if (wantr) {
					r[(n * k) + k] = 0.0f;
				}
			}
//!# Step 5: orthogonalization.
			if (k < n) {
				for (j = k + 1; j < n; j++) {
					dot = mc_dotpmx1f(m, n, n, k, j, q, q, 1);
					if (wantr) {
						r[(n * k) + j] = dot;
					}
					for (i = 0; i < m; i++) {
						q[(n * i) + j] = q[(n * i) + j] - (dot * q[(n * i) + k]);
					}
					if (wantpvi) {
						s    = mc_l2normmx1f(m, n, j, q);
						w[j] = mc_raise2f(s);
					}
				}
			} 
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_orthrmxnff(const int m, const int n, const float * a, float tol, double * q, double * MC_TARGET_RESTRICT r, double * MC_TARGET_RESTRICT w, int * pvi)
{
//!# Requires a[m x n], q[m x n] and r[n x n] if !null and w[n] && pvi[n] if !null where 1 < n <= m.
//!# Forming a ortho-normalized basis Q using Modified Gram-Schmidt method + a decimeting column step
//!# if norm < tol + iterative re-orthogonalization step for rank deficient systems. If R is not null
//!# upper-right-triangle is formed. @see Achiya Dax, `A modified Gram-schmidt algorithm with iterative
//!# orthogonalization and column pivoting`.
	const int wantr   = mc_nonnullptr(r);
	const int wantpvi = mc_nonnullptr(w) && mc_nonnullptr(pvi);

	int i, j, k, l, p;
	double bnorm, cnorm, dot, told, s;

	if (m >= n) {
		mc_copymxnff(m, n, q, a);

		if (wantr) {
			mc_eyenxn(n, r, 0);
		}

		if (wantpvi) {
			for (j = 0; j < n; j++) {
				s      = mc_l2normmx1(m, n, j, q);
				w[j]   = mc_raise2(s);
				pvi[j] = j;
			}
		}

		if (tol <= 0.0f) {
			tol = mc_fabsf(tol);
		}
		if (tol == 0.0f) {
			tol = MCLIMITS_EPSILONF;
		}
		told  = mc_cast(double, tol);
		bnorm = 0.0;
		for (k = 0; k < n; k++) {
//!# Step 1: pivoting if required.
			if (wantpvi) {
				mc_minmax1xn(n - k, w + k, MC_NULLPTR, &s, MC_NULLPTR, &l);
				l = l + k;
				if (k != l) {
					for (i = 0; i < m; i++) {
						mcswap_var(s, q[(n * i) + k], q[(n * i) + l]);
						if (i < k) {
							mcswap_var(s, r[(n * i) + k], r[(n * i) + l]);
						}
					}
					mcswap_var(p, pvi[k], pvi[l]);
				}
			}
//!# Step 2: re-orthogonalization.
			if (k > 0) {
				for (i = 0; i < k; i++) {
					dot = mc_dotpmx1(m, n, n, i, k, q, q, 1);
					if (wantr) {
						r[(n * i) + k] = r[(n * i) + k] + dot;
					}
					for (j = 0; j < m; j++) {
						q[(n * j) + k] = q[(n * j) + k] - (dot * q[(n * j) + i]);
					}
				}
			}
//!# Step 3: normalization.
			cnorm = mc_l2normmx1(m, n, k, q);
			if (cnorm != 0.0) {
//!# Step 4: close to zero decimation step.
				if (cnorm < told * bnorm) {
					mc_zerosmx1(m, n, k, q);
					q[(n * k) + k] = 1.0;
					if (wantr) {
						r[(n * k) + k] = 0.0;
					}
				} else {
					bnorm = mc_fmax(bnorm, cnorm);
					if (wantr) {
						r[(n * k) + k] = cnorm;
					}
					if (!mc_fisnear(cnorm, 1.0, 1)) {
						cnorm = 1.0 / cnorm;
						for (j = 0; j < m; j++) {
							q[(n * j) + k] = q[(n * j) + k] * cnorm;
						}
					}
				}
			} else {
				q[(n * k) + k] = 1.0;
				if (wantr) {
					r[(n * k) + k] = 0.0;
				}
			}
//!# Step 5: orthogonalization.
			if (k < n) {
				for (j = k + 1; j < n; j++) {
					dot = mc_dotpmx1(m, n, n, k, j, q, q, 1);
					if (wantr) {
						r[(n * k) + j] = dot;
					}
					for (i = 0; i < m; i++) {
						q[(n * i) + j] = q[(n * i) + j] - (dot * q[(n * i) + k]);
					}
					if (wantpvi) {
						s    = mc_l2normmx1(m, n, j, q);
						w[j] = mc_raise2(s);
					}
				}
			} 
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_orthrmxn(const int m, const int n, const double * a, double tol, double * q, double * MC_TARGET_RESTRICT r, double * MC_TARGET_RESTRICT w, int * pvi)
{
//!# Requires a[m x n], q[m x n] and r[n x n] if !null and w[n] && pvi[n] if !null where 1 < n <= m.
//!# A and Q may be the same. Forming a ortho-normalized basis Q using Modified Gram-Schmidt method
//!# + a decimeting column step if norm < tol + iterative re-orthogonalization step for rank deficient
//!# systems. If R is not null upper-right-triangle is formed. @see Achiya Dax, `A modified Gram-schmidt
//!# algorithm with iterative orthogonalization and column pivoting`.
	const int wantr   = mc_nonnullptr(r);
	const int wantpvi = mc_nonnullptr(w) && mc_nonnullptr(pvi);

	int i, j, k, l, p;
	double bnorm, cnorm, dot, s;

	if (m >= n) {
		if (a != q) {
			mc_copymxn(m, n, q, a);
		}
		if (wantr) {
			mc_eyenxn(n, r, 0);
		}

		if (wantpvi) {
			for (j = 0; j < n; j++) {
				s      = mc_l2normmx1(m, n, j, q);
				w[j]   = mc_raise2(s);
				pvi[j] = j;
			}
		}

		if (tol <= 0.0) {
			tol = mc_fabs(tol);
		}
		if (tol == 0.0) {
			tol = MCLIMITS_EPSILON;
		}
		bnorm = 0.0;
		for (k = 0; k < n; k++) {
//!# Step 1: pivoting if required.
			if (wantpvi) {
				mc_minmax1xn(n - k, w + k, MC_NULLPTR, &s, MC_NULLPTR, &l);
				l = l + k;
				if (k != l) {
					for (i = 0; i < m; i++) {
						mcswap_var(s, q[(n * i) + k], q[(n * i) + l]);
						if (i < k) {
							mcswap_var(s, r[(n * i) + k], r[(n * i) + l]);
						}
					}
					mcswap_var(p, pvi[k], pvi[l]);
				}
			}
//!# Step 2: re-orthogonalization.
			if (k > 0) {
				for (i = 0; i < k; i++) {
					dot = mc_dotpmx1(m, n, n, i, k, q, q, 1);
					if (wantr) {
						r[(n * i) + k] = r[(n * i) + k] + dot;
					}
					for (j = 0; j < m; j++) {
						q[(n * j) + k] = q[(n * j) + k] - (dot * q[(n * j) + i]);
					}
				}
			}
//!# Step 3: normalization.
			cnorm = mc_l2normmx1(m, n, k, q);
			if (cnorm != 0.0) {
//!# Step 4: close to zero decimation step.
				if (cnorm < tol * bnorm) {
					mc_zerosmx1(m, n, k, q);
					q[(n * k) + k] = 1.0;
					if (wantr) {
						r[(n * k) + k] = 0.0;
					}
				} else {
					bnorm = mc_fmax(bnorm, cnorm);
					if (wantr) {
						r[(n * k) + k] = cnorm;
					}
					if (!mc_fisnear(cnorm, 1.0, 1)) {
						cnorm = 1.0 / cnorm;
						for (j = 0; j < m; j++) {
							q[(n * j) + k] = q[(n * j) + k] * cnorm;
						}
					}
				}
			} else {
				q[(n * k) + k] = 1.0;
				if (wantr) {
					r[(n * k) + k] = 0.0;
				}
			}
//!# Step 5: orthogonalization.
			if (k < n) {
				for (j = k + 1; j < n; j++) {
					dot = mc_dotpmx1(m, n, n, k, j, q, q, 1);
					if (wantr) {
						r[(n * k) + j] = dot;
					}
					for (i = 0; i < m; i++) {
						q[(n * i) + j] = q[(n * i) + j] - (dot * q[(n * i) + k]);
					}
					if (wantpvi) {
						s    = mc_l2normmx1(m, n, j, q);
						w[j] = mc_raise2(s);
					}
				}
			} 
		}
		return 0;
	}
	return -1;
}

MC_TARGET_FUNC int mc_orthrmxnl(const int m, const int n, const long double * a, long double tol, long double * q, long double * MC_TARGET_RESTRICT r, long double * MC_TARGET_RESTRICT w, int * pvi)
{
//!# Requires a[m x n], q[m x n] and r[n x n] if !null and w[n] && pvi[n] if !null where 1 < n <= m.
//!# A and Q may be the same. Forming a ortho-normalized basis Q using Modified Gram-Schmidt method
//!# + a decimeting column step if norm < tol + iterative re-orthogonalization step for rank deficient
//!# systems. If R is not null upper-right-triangle is formed. @see Achiya Dax, `A modified Gram-schmidt
//!# algorithm with iterative orthogonalization and column pivoting`.
	const int wantr   = mc_nonnullptr(r);
	const int wantpvi = mc_nonnullptr(w) && mc_nonnullptr(pvi);

	int i, j, k, l, p;
	long double bnorm, cnorm, dot, s;

	if (m >= n) {
		if (a != q) {
			mc_copymxnl(m, n, q, a);
		}
		if (wantr) {
			mc_eyenxnl(n, r, 0);
		}

		if (wantpvi) {
			for (j = 0; j < n; j++) {
				s      = mc_l2normmx1l(m, n, j, q);
				w[j]   = mc_raise2l(s);
				pvi[j] = j;
			}
		}

		if (tol <= 0.0L) {
			tol = mc_fabsl(tol);
		}
		if (tol == 0.0) {
			tol = MCLIMITS_EPSILONL;
		}
		bnorm = 0.0L;
		for (k = 0; k < n; k++) {
//!# Step 1: pivoting if required.
			if (wantpvi) {
				mc_minmax1xnl(n - k, w + k, MC_NULLPTR, &s, MC_NULLPTR, &l);
				l = l + k;
				if (k != l) {
					for (i = 0; i < m; i++) {
						mcswap_var(s, q[(n * i) + k], q[(n * i) + l]);
						if (i < k) {
							mcswap_var(s, r[(n * i) + k], r[(n * i) + l]);
						}
					}
					mcswap_var(p, pvi[k], pvi[l]);
				}
			}
//!# Step 2: re-orthogonalization.
			if (k > 0) {
				for (i = 0; i < k; i++) {
					dot = mc_dotpmx1l(m, n, n, i, k, q, q, 1);
					if (wantr) {
						r[(n * i) + k] = r[(n * i) + k] + dot;
					}
					for (j = 0; j < m; j++) {
						q[(n * j) + k] = q[(n * j) + k] - (dot * q[(n * j) + i]);
					}
				}
			}
//!# Step 3: normalization.
			cnorm = mc_l2normmx1l(m, n, k, q);
			if (cnorm != 0.0L) {
//!# Step 4: close to zero decimation step.
				if (cnorm < tol * bnorm) {
					mc_zerosmx1l(m, n, k, q);
					q[(n * k) + k] = 1.0L;
					if (wantr) {
						r[(n * k) + k] = 0.0L;
					}
				} else {
					bnorm = mc_fmaxl(bnorm, cnorm);
					if (wantr) {
						r[(n * k) + k] = cnorm;
					}
					if (!mc_fisnearl(cnorm, 1.0L, 1)) {
						cnorm = 1.0L / cnorm;
						for (j = 0; j < m; j++) {
							q[(n * j) + k] = q[(n * j) + k] * cnorm;
						}
					}
				}
			} else {
				q[(n * k) + k] = 1.0L;
				if (wantr) {
					r[(n * k) + k] = 0.0L;
				}
			}
//!# Step 5: orthogonalization.
			if (k < n) {
				for (j = k + 1; j < n; j++) {
					dot = mc_dotpmx1l(m, n, n, k, j, q, q, 1);
					if (wantr) {
						r[(n * k) + j] = dot;
					}
					for (i = 0; i < m; i++) {
						q[(n * i) + j] = q[(n * i) + j] - (dot * q[(n * i) + k]);
					}
					if (wantpvi) {
						s    = mc_l2normmx1l(m, n, j, q);
						w[j] = mc_raise2l(s);
					}
				}
			} 
		}
		return 0;
	}
	return -1;
}

#endif /* !MC_ORTHRMXN_H */

/* EOF */
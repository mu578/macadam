//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lasq3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_lapack_lamch.h>
#include <macadam/lapack/mc_lapack_lasq4.h>
#include <macadam/lapack/mc_lapack_lasq5.h>
#include <macadam/lapack/mc_lapack_lasq6.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fmax.h>
#include <macadam/details/math/mc_fmin.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_LAPACKE_LASQ3_H
#define MC_LAPACKE_LASQ3_H

#pragma mark - mc_lapack_slasq3 -

MC_TARGET_PROC void mc_lapack_slasq3(const int i0, int * n0, float * z
	, int *    pp
	, float *  dmin
	, float *  sigma
	, float *  desig
	, float *  qmax
	, int *    nfail
	, int *    iter
	, int *     ndiv
	, const int ieee
	, int *     ttype
	, float *   dmin1
	, float *   dmin2
	, float *   dn
	, float *   dn1
	, float *   dn2
	, float *   g
	, float *   tau
) {
	const float hundrd = 100.0f, two = 2.0f, one = 1.0f, half = 0.5f, qurtr = 0.25f, zero = 0.0f;
	const float cbias = 1.5f;

	int ipn4, j4, n0in, nn;
	float eps, s, t, temp, tol, tol2;

	n0in = (*n0);
	eps  = mc_lapack_slamch('P');
	tol  = eps * hundrd;
	tol2 = mc_raise2f(tol);

F10:
	if ((*n0) < i0) {
		return;
	}
	if ((*n0) == i0) {
		goto F20;
	}
	nn = (4 * (*n0)) + (*pp);
	if ((*n0) == i0 + 1) {
		goto F40;
	}

	if (
		   mc_blas_vector_at(z, nn - 5) > tol2 * ((*sigma) + mc_blas_vector_at(z, nn - 3))
		&& mc_blas_vector_at(z, nn - (2 * (*pp)) - 4) > tol2 * mc_blas_vector_at(z, nn - 7)
	) {
		goto F30;
	}

F20:
	mc_blas_vector_at(z, (4 * (*n0)) - 3) = mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 3) + (*sigma);
	*n0                                   = (*n0) - 1;
	goto F10;

F30:
	if (
		   mc_blas_vector_at(z, nn - 9) > tol2 * (*sigma)
		&& mc_blas_vector_at(z, nn - (2 * (*pp)) - 8) > tol2 * mc_blas_vector_at(z, nn - 11)
	) {
		goto F50;
	}

F40:
	if (mc_blas_vector_at(z, nn - 3) > mc_blas_vector_at(z, nn - 7)) {
		s                            = mc_blas_vector_at(z, nn - 3);
		mc_blas_vector_at(z, nn - 3) = mc_blas_vector_at(z, nn - 7);
		mc_blas_vector_at(z, nn - 7) = s;
	}
	t = (mc_blas_vector_at(z, nn - 7) - mc_blas_vector_at(z, nn - 3) + mc_blas_vector_at(z, nn - 5)) * half;
	if (mc_blas_vector_at(z, nn - 5) > mc_blas_vector_at(z, nn - 3) * tol2 && t != zero) {
		s = mc_blas_vector_at(z, nn - 3) * (mc_blas_vector_at(z, nn - 5) / t);
		if (s <= t) {
			s = mc_blas_vector_at(z, nn - 3) * (mc_blas_vector_at(z, nn - 5) / (t * (mc_sqrtf(s / t + one) + one)));
		} else {
			s = mc_blas_vector_at(z, nn - 3) * (mc_blas_vector_at(z, nn - 5) / (t + mc_sqrtf(t) * mc_sqrtf(t + s)));
		}
		t                            = mc_blas_vector_at(z, nn - 7) + (s + mc_blas_vector_at(z, nn - 5));
		mc_blas_vector_at(z, nn - 3) = mc_blas_vector_at(z, nn - 3) * (mc_blas_vector_at(z, nn - 7) / t);
		mc_blas_vector_at(z, nn - 7) = t;
	}
	mc_blas_vector_at(z, (4 * (*n0)) - 7) = mc_blas_vector_at(z, nn - 7) + (*sigma);
	mc_blas_vector_at(z, (4 * (*n0)) - 3) = mc_blas_vector_at(z, nn - 3) + (*sigma);
	*n0                                   = (*n0) - 2;
	goto F10;

F50:
	if ((*pp)== 2) {
		*pp = 0;
	}

	if (*dmin <= zero || *n0 < n0in) {
		if (mc_blas_vector_at(z, (4 * i0) + (*pp) - 3) * cbias < mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 3)) {
			 ipn4                                        = 4 * (i0 + (*n0));
			for (j4 = (4 * i0); j4 <= (2 * (i0 + (*n0) - 1 )); j4 += 4) {
				temp                                = mc_blas_vector_at(z, j4 - 3);
				mc_blas_vector_at(z, j4 - 3)        = mc_blas_vector_at(z, ipn4 - j4 - 3);
				mc_blas_vector_at(z, ipn4 - j4 - 3) = temp;
				temp                                = mc_blas_vector_at(z, j4 - 2);
				mc_blas_vector_at(z, j4 - 2)        = mc_blas_vector_at(z, ipn4 - j4 - 2);
				mc_blas_vector_at(z, ipn4 - j4 - 2) = temp;
				temp                                = mc_blas_vector_at(z, j4 - 1);
				mc_blas_vector_at(z, j4 - 1)        = mc_blas_vector_at(z, ipn4 - j4 - 5);
				mc_blas_vector_at(z, ipn4 - j4 - 5) = temp;
				temp                                = mc_blas_vector_at(z, j4);
				mc_blas_vector_at(z, j4)            = mc_blas_vector_at(z, ipn4 - j4 - 4);
				mc_blas_vector_at(z, ipn4 - j4 - 4) = temp;
			}
			if (*n0 - i0 <= 4) {
				mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 1) = mc_blas_vector_at(z, (4 * i0) + (*pp) - 1);
				mc_blas_vector_at(z, (4 * (*n0)) - (*pp))     = mc_blas_vector_at(z, (4 * i0) - *pp);
			}
			*dmin2                                        = mc_fminf(*dmin2, mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 1));
			mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 1) = mc_fminf(
				  mc_fminf(mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 1), mc_blas_vector_at(z, (4 * i0) + (*pp) - 1))
				, mc_blas_vector_at(z, (4 * i0) + (*pp) + 3)
			);
			mc_blas_vector_at(z, (4 * (*n0)) - (*pp))     = mc_fminf(
				  mc_fminf(mc_blas_vector_at(z, (4 * (*n0)) - (*pp)) , mc_blas_vector_at(z, (4 * i0) - *pp))
				, mc_blas_vector_at(z, (4 * i0) - (*pp) + 4)
			);
			*qmax                                         = mc_fmaxf(
				  mc_fmaxf(*qmax, mc_blas_vector_at(z, (4 * i0) + (*pp) - 3))
				, mc_blas_vector_at(z, (4 * i0) + (*pp) + 1)
			);
			*dmin                                         = -zero;
		}
	}

	mc_lapack_slasq4(i0, *n0, z, *pp, n0in, *dmin, *dmin1, *dmin2, *dn, *dn1, *dn2, tau, ttype, g);

F70:
	mc_lapack_slasq5(i0, *n0, z, *pp, *tau, *sigma, dmin, dmin1, dmin2, dn, dn1, dn2, ieee, eps);
	*ndiv = (*ndiv) + (*n0 - i0 + 2);
	*iter = (*iter) + 1;

	if (*dmin >= zero && *dmin1 >= zero) {
		goto F90;
	} else if (
		   *dmin < zero && *dmin1 > zero
		&& mc_blas_vector_at(z, (4 * (*n0 - 1)) - (*pp)) < tol * ((*sigma) + (*dn1))
		&& mc_fabsf(*dn) < tol * (*sigma)
	) {
		mc_blas_vector_at(z, (4 * (*n0 - 1) - (*pp) + 2)) = zero;
		*dmin                                             = zero;
		goto F90;
	} else if (*dmin < zero) {
		*nfail = (*nfail) + 1;
		if (*ttype < -22) {
			*tau = zero;
		} else if (*dmin1 > zero) {
			*tau   = (*tau + *dmin) * (one - eps * two);
			*ttype = *ttype + -11;
		} else {
			*tau   = *tau * qurtr;
			*ttype = *ttype + -12;
		}
		goto F70;
	} else if (mc_isnan(*dmin)) {
		if (*tau == zero) {
			goto F80;
		} else {
			*tau = zero;
			goto F70;
		}
	} else {
		goto F80;
	}

F80:
	mc_lapack_slasq6(i0, *n0, z, *pp, dmin, dmin1, dmin2, dn, dn1, dn2);
	*ndiv = (*ndiv) + ((*n0) - i0 + 2);
	*iter = (*iter) + 1;
	*tau  = zero;

F90:
	if (*tau < (*sigma)) {
		*desig = (*desig) + (*tau);
		 t     = (*sigma) + (*desig);
		*desig = (*desig) - (t - (*sigma));
	} else {
		 t     = (*sigma) + (*tau);
		*desig = (*sigma) - (t - (*tau)) + (*desig);
	}
	*sigma = t;
}

#pragma mark - mc_lapack_dlasq3 -

MC_TARGET_PROC void mc_lapack_dlasq3(const int i0, int * n0, double * z
	, int *     pp
	, double *  dmin
	, double *  sigma
	, double *  desig
	, double *  qmax
	, int *     nfail
	, int *     iter
	, int *     ndiv
	, const int ieee
	, int *     ttype
	, double *  dmin1
	, double *  dmin2
	, double *  dn
	, double *  dn1
	, double *  dn2
	, double *  g
	, double *  tau
) {
	const double hundrd = 100.0, two = 2.0, one = 1.0, half = 0.5, qurtr = 0.25, zero = 0.0;
	const double cbias = 1.5;

	int ipn4, j4, n0in, nn;
	double eps, s, t, temp, tol, tol2;

	n0in = (*n0);
	eps  = mc_lapack_dlamch('P');
	tol  = eps * hundrd;
	tol2 = mc_raise2(tol);

F10:
	if ((*n0) < i0) {
		return;
	}
	if ((*n0) == i0) {
		goto F20;
	}
	nn = (4 * (*n0)) + (*pp);
	if ((*n0) == i0 + 1) {
		goto F40;
	}

	if (
		   mc_blas_vector_at(z, nn - 5) > tol2 * ((*sigma) + mc_blas_vector_at(z, nn - 3))
		&& mc_blas_vector_at(z, nn - (2 * (*pp)) - 4) > tol2 * mc_blas_vector_at(z, nn - 7)
	) {
		goto F30;
	}

F20:
	mc_blas_vector_at(z, (4 * (*n0)) - 3) = mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 3) + (*sigma);
	*n0                                   = (*n0) - 1;
	goto F10;

F30:
	if (
		   mc_blas_vector_at(z, nn - 9) > tol2 * (*sigma)
		&& mc_blas_vector_at(z, nn - (2 * (*pp)) - 8) > tol2 * mc_blas_vector_at(z, nn - 11)
	) {
		goto F50;
	}

F40:
	if (mc_blas_vector_at(z, nn - 3) > mc_blas_vector_at(z, nn - 7)) {
		s                            = mc_blas_vector_at(z, nn - 3);
		mc_blas_vector_at(z, nn - 3) = mc_blas_vector_at(z, nn - 7);
		mc_blas_vector_at(z, nn - 7) = s;
	}
	t = (mc_blas_vector_at(z, nn - 7) - mc_blas_vector_at(z, nn - 3) + mc_blas_vector_at(z, nn - 5)) * half;
	if (mc_blas_vector_at(z, nn - 5) > mc_blas_vector_at(z, nn - 3) * tol2 && t != zero) {
		s = mc_blas_vector_at(z, nn - 3) * (mc_blas_vector_at(z, nn - 5) / t);
		if (s <= t) {
			s = mc_blas_vector_at(z, nn - 3) * (mc_blas_vector_at(z, nn - 5) / (t * (mc_sqrt(s / t + one) + one)));
		} else {
			s = mc_blas_vector_at(z, nn - 3) * (mc_blas_vector_at(z, nn - 5) / (t + mc_sqrt(t) * mc_sqrt(t + s)));
		}
		t                            = mc_blas_vector_at(z, nn - 7) + (s + mc_blas_vector_at(z, nn - 5));
		mc_blas_vector_at(z, nn - 3) = mc_blas_vector_at(z, nn - 3) * (mc_blas_vector_at(z, nn - 7) / t);
		mc_blas_vector_at(z, nn - 7) = t;
	}
	mc_blas_vector_at(z, (4 * (*n0)) - 7) = mc_blas_vector_at(z, nn - 7) + (*sigma);
	mc_blas_vector_at(z, (4 * (*n0)) - 3) = mc_blas_vector_at(z, nn - 3) + (*sigma);
	*n0                                   = (*n0) - 2;
	goto F10;

F50:
	if ((*pp)== 2) {
		*pp = 0;
	}

	if (*dmin <= zero || *n0 < n0in) {
		if (mc_blas_vector_at(z, (4 * i0) + (*pp) - 3) * cbias < mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 3)) {
			 ipn4                                        = 4 * (i0 + (*n0));
			for (j4 = (4 * i0); j4 <= (2 * (i0 + (*n0) - 1 )); j4 += 4) {
				temp                                = mc_blas_vector_at(z, j4 - 3);
				mc_blas_vector_at(z, j4 - 3)        = mc_blas_vector_at(z, ipn4 - j4 - 3);
				mc_blas_vector_at(z, ipn4 - j4 - 3) = temp;
				temp                                = mc_blas_vector_at(z, j4 - 2);
				mc_blas_vector_at(z, j4 - 2)        = mc_blas_vector_at(z, ipn4 - j4 - 2);
				mc_blas_vector_at(z, ipn4 - j4 - 2) = temp;
				temp                                = mc_blas_vector_at(z, j4 - 1);
				mc_blas_vector_at(z, j4 - 1)        = mc_blas_vector_at(z, ipn4 - j4 - 5);
				mc_blas_vector_at(z, ipn4 - j4 - 5) = temp;
				temp                                = mc_blas_vector_at(z, j4);
				mc_blas_vector_at(z, j4)            = mc_blas_vector_at(z, ipn4 - j4 - 4);
				mc_blas_vector_at(z, ipn4 - j4 - 4) = temp;
			}
			if (*n0 - i0 <= 4) {
				mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 1) = mc_blas_vector_at(z, (4 * i0) + (*pp) - 1);
				mc_blas_vector_at(z, (4 * (*n0)) - (*pp))     = mc_blas_vector_at(z, (4 * i0) - *pp);
			}
			*dmin2                                        = mc_fmin(*dmin2, mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 1));
			mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 1) = mc_fmin(
				  mc_fmin(mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 1), mc_blas_vector_at(z, (4 * i0) + (*pp) - 1))
				, mc_blas_vector_at(z, (4 * i0) + (*pp) + 3)
			);
			mc_blas_vector_at(z, (4 * (*n0)) - (*pp))     = mc_fmin(
				  mc_fmin(mc_blas_vector_at(z, (4 * (*n0)) - (*pp)) , mc_blas_vector_at(z, (4 * i0) - *pp))
				, mc_blas_vector_at(z, (4 * i0) - (*pp) + 4)
			);
			*qmax                                         = mc_fmax(
				  mc_fmax(*qmax, mc_blas_vector_at(z, (4 * i0) + (*pp) - 3))
				, mc_blas_vector_at(z, (4 * i0) + (*pp) + 1)
			);
			*dmin                                         = -zero;
		}
	}

	mc_lapack_dlasq4(i0, *n0, z, *pp, n0in, *dmin, *dmin1, *dmin2, *dn, *dn1, *dn2, tau, ttype, g);

F70:
	mc_lapack_dlasq5(i0, *n0, z, *pp, *tau, *sigma, dmin, dmin1, dmin2, dn, dn1, dn2, ieee, eps);
	*ndiv = (*ndiv) + (*n0 - i0 + 2);
	*iter = (*iter) + 1;

	if (*dmin >= zero && *dmin1 >= zero) {
		goto F90;
	} else if (
		   *dmin < zero && *dmin1 > zero
		&& mc_blas_vector_at(z, (4 * (*n0 - 1)) - (*pp)) < tol * ((*sigma) + (*dn1))
		&& mc_fabs(*dn) < tol * (*sigma)
	) {
		mc_blas_vector_at(z, (4 * (*n0 - 1) - (*pp) + 2)) = zero;
		*dmin                                             = zero;
		goto F90;
	} else if (*dmin < zero) {
		*nfail = (*nfail) + 1;
		if (*ttype < -22) {
			*tau = zero;
		} else if (*dmin1 > zero) {
			*tau   = (*tau + *dmin) * (one - eps * two);
			*ttype = *ttype + -11;
		} else {
			*tau   = *tau * qurtr;
			*ttype = *ttype + -12;
		}
		goto F70;
	} else if (mc_isnan(*dmin)) {
		if (*tau == zero) {
			goto F80;
		} else {
			*tau = zero;
			goto F70;
		}
	} else {
		goto F80;
	}

F80:
	mc_lapack_dlasq6(i0, *n0, z, *pp, dmin, dmin1, dmin2, dn, dn1, dn2);
	*ndiv = (*ndiv) + ((*n0) - i0 + 2);
	*iter = (*iter) + 1;
	*tau  = zero;

F90:
	if (*tau < (*sigma)) {
		*desig = (*desig) + (*tau);
		 t     = (*sigma) + (*desig);
		*desig = (*desig) - (t - (*sigma));
	} else {
		 t     = (*sigma) + (*tau);
		*desig = (*sigma) - (t - (*tau)) + (*desig);
	}
	*sigma = t;
}

#pragma mark - mc_lapack_llasq3 -

MC_TARGET_PROC void mc_lapack_llasq3(const int i0, int * n0, long double * z
	, int *         pp
	, long double * dmin
	, long double * sigma
	, long double * desig
	, long double * qmax
	, int *         nfail
	, int *         iter
	, int *         ndiv
	, const int     ieee
	, int *         ttype
	, long double * dmin1
	, long double * dmin2
	, long double * dn
	, long double * dn1
	, long double * dn2
	, long double * g
	, long double * tau
) {
	const long double hundrd = 100.0L, two = 2.0L, one = 1.0L, half = 0.5L, qurtr = 0.25L, zero = 0.0L;
	const long double cbias = 1.5L;

	int ipn4, j4, n0in, nn;
	long double eps, s, t, temp, tol, tol2;

	n0in = (*n0);
	eps  = mc_lapack_llamch('P');
	tol  = eps * hundrd;
	tol2 = mc_raise2l(tol);

F10:
	if ((*n0) < i0) {
		return;
	}
	if ((*n0) == i0) {
		goto F20;
	}
	nn = (4 * (*n0)) + (*pp);
	if ((*n0) == i0 + 1) {
		goto F40;
	}

	if (
		   mc_blas_vector_at(z, nn - 5) > tol2 * ((*sigma) + mc_blas_vector_at(z, nn - 3))
		&& mc_blas_vector_at(z, nn - (2 * (*pp)) - 4) > tol2 * mc_blas_vector_at(z, nn - 7)
	) {
		goto F30;
	}

F20:
	mc_blas_vector_at(z, (4 * (*n0)) - 3) = mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 3) + (*sigma);
	*n0                                   = (*n0) - 1;
	goto F10;

F30:
	if (
		   mc_blas_vector_at(z, nn - 9) > tol2 * (*sigma)
		&& mc_blas_vector_at(z, nn - (2 * (*pp)) - 8) > tol2 * mc_blas_vector_at(z, nn - 11)
	) {
		goto F50;
	}

F40:
	if (mc_blas_vector_at(z, nn - 3) > mc_blas_vector_at(z, nn - 7)) {
		s                            = mc_blas_vector_at(z, nn - 3);
		mc_blas_vector_at(z, nn - 3) = mc_blas_vector_at(z, nn - 7);
		mc_blas_vector_at(z, nn - 7) = s;
	}
	t = (mc_blas_vector_at(z, nn - 7) - mc_blas_vector_at(z, nn - 3) + mc_blas_vector_at(z, nn - 5)) * half;
	if (mc_blas_vector_at(z, nn - 5) > mc_blas_vector_at(z, nn - 3) * tol2 && t != zero) {
		s = mc_blas_vector_at(z, nn - 3) * (mc_blas_vector_at(z, nn - 5) / t);
		if (s <= t) {
			s = mc_blas_vector_at(z, nn - 3) * (mc_blas_vector_at(z, nn - 5) / (t * (mc_sqrtl(s / t + one) + one)));
		} else {
			s = mc_blas_vector_at(z, nn - 3) * (mc_blas_vector_at(z, nn - 5) / (t + mc_sqrtl(t) * mc_sqrtl(t + s)));
		}
		t                            = mc_blas_vector_at(z, nn - 7) + (s + mc_blas_vector_at(z, nn - 5));
		mc_blas_vector_at(z, nn - 3) = mc_blas_vector_at(z, nn - 3) * (mc_blas_vector_at(z, nn - 7) / t);
		mc_blas_vector_at(z, nn - 7) = t;
	}
	mc_blas_vector_at(z, (4 * (*n0)) - 7) = mc_blas_vector_at(z, nn - 7) + (*sigma);
	mc_blas_vector_at(z, (4 * (*n0)) - 3) = mc_blas_vector_at(z, nn - 3) + (*sigma);
	*n0                                   = (*n0) - 2;
	goto F10;

F50:
	if ((*pp)== 2) {
		*pp = 0;
	}

	if (*dmin <= zero || *n0 < n0in) {
		if (mc_blas_vector_at(z, (4 * i0) + (*pp) - 3) * cbias < mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 3)) {
			 ipn4                                        = 4 * (i0 + (*n0));
			for (j4 = (4 * i0); j4 <= (2 * (i0 + (*n0) - 1 )); j4 += 4) {
				temp                                = mc_blas_vector_at(z, j4 - 3);
				mc_blas_vector_at(z, j4 - 3)        = mc_blas_vector_at(z, ipn4 - j4 - 3);
				mc_blas_vector_at(z, ipn4 - j4 - 3) = temp;
				temp                                = mc_blas_vector_at(z, j4 - 2);
				mc_blas_vector_at(z, j4 - 2)        = mc_blas_vector_at(z, ipn4 - j4 - 2);
				mc_blas_vector_at(z, ipn4 - j4 - 2) = temp;
				temp                                = mc_blas_vector_at(z, j4 - 1);
				mc_blas_vector_at(z, j4 - 1)        = mc_blas_vector_at(z, ipn4 - j4 - 5);
				mc_blas_vector_at(z, ipn4 - j4 - 5) = temp;
				temp                                = mc_blas_vector_at(z, j4);
				mc_blas_vector_at(z, j4)            = mc_blas_vector_at(z, ipn4 - j4 - 4);
				mc_blas_vector_at(z, ipn4 - j4 - 4) = temp;
			}
			if (*n0 - i0 <= 4) {
				mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 1) = mc_blas_vector_at(z, (4 * i0) + (*pp) - 1);
				mc_blas_vector_at(z, (4 * (*n0)) - (*pp))     = mc_blas_vector_at(z, (4 * i0) - *pp);
			}
			*dmin2                                        = mc_fminl(*dmin2, mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 1));
			mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 1) = mc_fminl(
				  mc_fminl(mc_blas_vector_at(z, (4 * (*n0)) + (*pp) - 1), mc_blas_vector_at(z, (4 * i0) + (*pp) - 1))
				, mc_blas_vector_at(z, (4 * i0) + (*pp) + 3)
			);
			mc_blas_vector_at(z, (4 * (*n0)) - (*pp))     = mc_fminl(
				  mc_fminl(mc_blas_vector_at(z, (4 * (*n0)) - (*pp)) , mc_blas_vector_at(z, (4 * i0) - *pp))
				, mc_blas_vector_at(z, (4 * i0) - (*pp) + 4)
			);
			*qmax                                         = mc_fmaxl(
				  mc_fmaxl(*qmax, mc_blas_vector_at(z, (4 * i0) + (*pp) - 3))
				, mc_blas_vector_at(z, (4 * i0) + (*pp) + 1)
			);
			*dmin                                         = -zero;
		}
	}

	mc_lapack_llasq4(i0, *n0, z, *pp, n0in, *dmin, *dmin1, *dmin2, *dn, *dn1, *dn2, tau, ttype, g);

F70:
	mc_lapack_llasq5(i0, *n0, z, *pp, *tau, *sigma, dmin, dmin1, dmin2, dn, dn1, dn2, ieee, eps);
	*ndiv = (*ndiv) + (*n0 - i0 + 2);
	*iter = (*iter) + 1;

	if (*dmin >= zero && *dmin1 >= zero) {
		goto F90;
	} else if (
		   *dmin < zero && *dmin1 > zero
		&& mc_blas_vector_at(z, (4 * (*n0 - 1)) - (*pp)) < tol * ((*sigma) + (*dn1))
		&& mc_fabsl(*dn) < tol * (*sigma)
	) {
		mc_blas_vector_at(z, (4 * (*n0 - 1) - (*pp) + 2)) = zero;
		*dmin                                             = zero;
		goto F90;
	} else if (*dmin < zero) {
		*nfail = (*nfail) + 1;
		if (*ttype < -22) {
			*tau = zero;
		} else if (*dmin1 > zero) {
			*tau   = (*tau + *dmin) * (one - eps * two);
			*ttype = *ttype + -11;
		} else {
			*tau   = *tau * qurtr;
			*ttype = *ttype + -12;
		}
		goto F70;
	} else if (mc_isnan(*dmin)) {
		if (*tau == zero) {
			goto F80;
		} else {
			*tau = zero;
			goto F70;
		}
	} else {
		goto F80;
	}

F80:
	mc_lapack_llasq6(i0, *n0, z, *pp, dmin, dmin1, dmin2, dn, dn1, dn2);
	*ndiv = (*ndiv) + ((*n0) - i0 + 2);
	*iter = (*iter) + 1;
	*tau  = zero;

F90:
	if (*tau < (*sigma)) {
		*desig = (*desig) + (*tau);
		 t     = (*sigma) + (*desig);
		*desig = (*desig) - (t - (*sigma));
	} else {
		 t     = (*sigma) + (*tau);
		*desig = (*sigma) - (t - (*tau)) + (*desig);
	}
	*sigma = t;
}

#endif /* !MC_LAPACKE_LASQ3_H */

/* EOF */
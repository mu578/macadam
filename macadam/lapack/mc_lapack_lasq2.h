//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lasq2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_lapack_lamch.h>
#include <macadam/lapack/mc_lapack_lasrt.h>
#include <macadam/lapack/mc_lapack_lasq3.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fmax.h>
#include <macadam/details/math/mc_fmin.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_LAPACKE_LASQ2_H
#define MC_LAPACKE_LASQ2_H

#pragma mark - mc_lapack_slasq2 -

MC_TARGET_PROC void mc_lapack_slasq2(const int n, float * z, int * info)
{
	const float hundrd = 100.0f, four = 4.0f, two = 2.0f, one = 1.0f, half = 0.5f, zero = 0.0f;
	const float cbias = 1.5f;

	int ieee;
	int i0, i4, iinfo, ipn4, iter, iwhila, iwhilb, k,
	    kmin, n0, nbig, ndiv, nfail, pp, splt, ttype,
	    i1, n1;
	float d, dee, deemin, desig, dmin, dmin1, dmin2, dn,
	      dn1, dn2, e, emax, emin, eps, g, oldemn, qmax,
	      qmin, s, safmin, sigma, t, tau, temp, tol,
	      tol2, trace, zmax, tempe, tempq;

	*info   = 0;
	 eps    = mc_lapack_slamch('P');
	 safmin = mc_lapack_slamch('S');
	 tol    = eps * hundrd;
	 tol2   = mc_raise2f(tol);

	if (n < 0) {
		*info = -1;
		mc_blas_xerbla("SLASQ2", 1);
		return;
	} else if (n == 0) {
		return;
	} else if (n == 1) {
		if (mc_blas_vector_at(z, 1) < zero) {
			*info = -201;
			mc_blas_xerbla("SLASQ2", 2);
		}
		return;
	} else if (n == 2) {
		if (mc_blas_vector_at(z, 2) < zero || mc_blas_vector_at(z, 3) < zero) {
			*info = -2;
			mc_blas_xerbla("SLASQ2", 2);
			return;
		} else if (mc_blas_vector_at(z, 3) > mc_blas_vector_at(z, 1)) {
			d                       = mc_blas_vector_at(z, 3);
			mc_blas_vector_at(z, 3) = mc_blas_vector_at(z, 1);
			mc_blas_vector_at(z, 1) = d;
		}
		mc_blas_vector_at(z, 5) = mc_blas_vector_at(z, 1) + mc_blas_vector_at(z, 2) + mc_blas_vector_at(z, 3);
		if (mc_blas_vector_at(z, 2) > mc_blas_vector_at(z, 3) * tol2) {
			t = (mc_blas_vector_at(z, 1) - mc_blas_vector_at(z, 3) + mc_blas_vector_at(z, 2)) * half;
			s = mc_blas_vector_at(z, 3) * (mc_blas_vector_at(z, 2) / t);
			if (s <= t) {
				s = mc_blas_vector_at(z, 3) * (mc_blas_vector_at(z, 2) / (t * (mc_sqrtf(s / t + one) + one)));
			} else {
				s = mc_blas_vector_at(z, 3) * (mc_blas_vector_at(z, 2) / (t + mc_sqrtf(t) * mc_sqrtf(t + s)));
			}
			t                       = mc_blas_vector_at(z, 1) + (s + mc_blas_vector_at(z, 2));
			mc_blas_vector_at(z, 3) = mc_blas_vector_at(z, 3) * (mc_blas_vector_at(z, 1) / t);
			mc_blas_vector_at(z, 1) = t;
		}
		mc_blas_vector_at(z, 2) = mc_blas_vector_at(z, 3);
		mc_blas_vector_at(z, 6) = mc_blas_vector_at(z, 2) + mc_blas_vector_at(z, 1);
		return;
	}

	mc_blas_vector_at(z, 2 * n) = zero;
	emin                        = mc_blas_vector_at(z, 2);
	qmax                        = zero;
	zmax                        = zero;
	d                           = zero;
	e                           = zero;

	for (k = 1; k <= (2 * (n - 1)); k += 2) {
		if (mc_blas_vector_at(z, k) < zero) {
			*info = -(200 + k);
			mc_blas_xerbla("SLASQ2", 2);
			return;
		} else if (mc_blas_vector_at(z, k + 1) < zero) {
			*info = -(200 + k + 1);
			mc_blas_xerbla("SLASQ2", 2);
			return;
		}
		d    = d + mc_blas_vector_at(z, k);
		e    = e + mc_blas_vector_at(z, k + 1);
		qmax = mc_fmaxf(qmax, mc_blas_vector_at(z, k));
		emin = mc_fminf(emin, mc_blas_vector_at(z, k + 1));
		zmax = mc_fmaxf(mc_maxmag(qmax, zmax), mc_blas_vector_at(z, k + 1));
	}
	if (mc_blas_vector_at(z, (2 * n) - 1) < zero) {
		*info = -(200 + ((2 * n) - 1));
		mc_blas_xerbla("SLASQ2", 2);
		return;
	}
	d    = d + mc_blas_vector_at(z, (2 * n) - 1);
	qmax = mc_fmaxf(qmax, mc_blas_vector_at(z, (2 * n) - 1));
	zmax = mc_fmaxf(qmax, zmax);
	mc_unused(zmax);

	if (e == zero) {
		for (k = 2; k <= n; ++k) {
			mc_blas_vector_at(z, k) = mc_blas_vector_at(z, (2 * k) - 1);
		}
		mc_lapack_slasrt('D', n, z, &iinfo);
		mc_blas_vector_at(z, (2 * n) - 1) = d;
		return;
	}

	trace = d + e;

	if (trace == zero) {
		mc_blas_vector_at(z, (2 * n) - 1) = zero;
		return;
	}

	ieee = 0;

	for (k = (2 * n); k >= 2; k += -2) {
		mc_blas_vector_at(z,  2 * k     ) = zero;
		mc_blas_vector_at(z, (2 * k) - 1) = mc_blas_vector_at(z, k);
		mc_blas_vector_at(z, (2 * k) - 2) = zero;
		mc_blas_vector_at(z, (2 * k) - 3) = mc_blas_vector_at(z, k - 1);
	}

	i0 = 1;
	n0 = n;

	if (mc_blas_vector_at(z, (4 * i0) - 3) * cbias < mc_blas_vector_at(z, (4 * n0) - 3)) {
		ipn4 = 4 * (i0 + n0);
		for (i4 = (4 * i0); i4 <= (2 * (i0 + n0 - 1)); i4 += 4) {
			temp                                = mc_blas_vector_at(z, i4 - 3);
			mc_blas_vector_at(z, i4 - 3)        = mc_blas_vector_at(z, ipn4 - i4 - 3);
			mc_blas_vector_at(z, ipn4 - i4 - 3) = temp;
			temp                                = mc_blas_vector_at(z, i4 - 1);
			mc_blas_vector_at(z, i4 - 1)        = mc_blas_vector_at(z, ipn4 - i4 - 5);
			mc_blas_vector_at(z, ipn4 - i4 - 5) = temp;
		}
	}

	pp = 0;

	for (k = 1; k <= 2; ++k) {
		d = mc_blas_vector_at(z, (4 * n0) + pp - 3);
		for (i4 = ((4 * (n0 - 1)) + pp); i4 >=  ((4 * i0) + pp); i4 += -4) {
			if (mc_blas_vector_at(z, i4 - 1) <= tol2 * d) {
				mc_blas_vector_at(z, i4 - 1) = -zero;
				d                            = mc_blas_vector_at(z, i4 - 3);
			} else {
				d = mc_blas_vector_at(z, i4 - 3) * (d / (d + mc_blas_vector_at(z, i4 - 1)));
			}
		}

		emin = mc_blas_vector_at(z, (4 * i0) + pp + 1);
		d    = mc_blas_vector_at(z, (4 * i0) + pp - 3);
		for (i4 = ((4 * i0) + pp); i4 <= ((4 * (n0 - 1)) + pp); i4 += 4) {
			mc_blas_vector_at(z, i4 - (2 * pp) - 2) = d + mc_blas_vector_at(z, i4 - 1);
			if (mc_blas_vector_at(z, i4 - 1) <= tol2 * d) {
				mc_blas_vector_at(z, i4 - 1)            = -zero;
				mc_blas_vector_at(z, i4 - (2 * pp) - 2) = d;
				mc_blas_vector_at(z, i4 - (2 * pp))     = zero;
				d                                       = mc_blas_vector_at(z, i4 + 1);
			} else if (
					safmin * mc_blas_vector_at(z, i4 + 1) < mc_blas_vector_at(z, i4 - (2 * pp) - 2)
				&& safmin * mc_blas_vector_at(z, i4 - (2 * pp) - 2) < mc_blas_vector_at(z, i4 + 1)
			) {
				temp                                = mc_blas_vector_at(z, i4 + 1) / mc_blas_vector_at(z, i4 - (2 * pp) - 2);
				mc_blas_vector_at(z, i4 - (2 * pp)) = mc_blas_vector_at(z, i4 - 1) * temp;
				d                                   = d * temp;
			} else {
				mc_blas_vector_at(z, i4 - (2 * pp)) = mc_blas_vector_at(z, i4 + 1) * (mc_blas_vector_at(z, i4 - 1) / mc_blas_vector_at(z, i4 - (2 * pp) - 2));
				d                                   = mc_blas_vector_at(z, i4 + 1) * (d / mc_blas_vector_at(z, i4 - (2 * pp) - 2));
			}
			emin = mc_fminf(emin, mc_blas_vector_at(z, i4 - (2 * pp)));
		}
		mc_blas_vector_at(z, (4 * n0) - pp - 2) = d;

		qmax = mc_blas_vector_at(z, (4 * i0) - pp - 2);
		for (i4 = ((4 * i0) - pp + 2); i4 <= ((4 * n0) - pp - 2); i4 += 4) {
			qmax = mc_fmaxf(qmax, mc_blas_vector_at(z, i4));
		}
		pp = 1 - pp;
	}

	ttype = 0;
	dmin1 = zero;
	dmin2 = zero;
	dn    = zero;
	dn1   = zero;
	dn2   = zero;
	g     = zero;
	tau   = zero;

	iter  = 2;
	nfail = 0;
	ndiv  = 2 * (n0 - i0);

	for (iwhila = 1; iwhila <= (n + 1); ++iwhila) {
		if (n0 < 1) {
			goto F170;
		}
		desig = zero;
		if (n0 == n) {
			sigma = zero;
		} else {
			sigma = -mc_blas_vector_at(z, (4 * n0) - 1);
		}
		if (sigma < zero) {
			*info = 1;
			return;
		}
		emax = zero;
		if (n0 > i0) {
			emin = mc_fabsf(mc_blas_vector_at(z, (4 * n0) - 5));
		} else {
			emin = zero;
		}
		qmin = mc_blas_vector_at(z, (4 * n0) - 3);
		qmax = qmin;
		for (i4 = (4 * n0); i4 >= 8; i4 += -4) {
			if (mc_blas_vector_at(z, i4 - 5) <= zero) {
				goto F100;
			}
			if (qmin >= emax * four) {
				qmin = mc_fminf(qmin, mc_blas_vector_at(z, i4 - 3));
				emax = mc_fmaxf(emax, mc_blas_vector_at(z, i4 - 5));
			}
			qmax = mc_fmaxf(qmax, mc_blas_vector_at(z, i4 - 7) + mc_blas_vector_at(z, i4 - 5));
			emin = mc_fminf(emin, mc_blas_vector_at(z, i4 - 5));
		}
		i4 = 4;

F100:
		i0 = i4 / 4;
		pp = 0;

		if ((n0 - i0) > 1) {
			dee    = mc_blas_vector_at(z, (4 * i0) - 3);
			deemin = dee;
			kmin   = i0;
			for (i4 = ((4 * i0) + 1); i4 <= ((4 * n0) - 3); i4 += 4) {
				dee = mc_blas_vector_at(z, i4) * (dee / (dee + mc_blas_vector_at(z, i4 - 2)));
				if (dee <= deemin) {
					deemin = dee;
					kmin = (i4 + 3) / 4;
				}
			}
			if (((kmin - i0) * 2) < n0 - kmin && deemin <= mc_blas_vector_at(z, (4 * n0) - 3) * half) {
				ipn4 = 4 * (i0 + n0);
				pp   = 2;
				for (i4 = (4 * i0); i4 <= (2 * (i0 + n0 - 1)); i4 += 4) {
					temp                                = mc_blas_vector_at(z, i4 - 3);
					mc_blas_vector_at(z, i4 - 3)        = mc_blas_vector_at(z, ipn4 - i4 - 3);
					mc_blas_vector_at(z, ipn4 - i4 - 3) = temp;
					temp                                = mc_blas_vector_at(z, i4 - 2);
					mc_blas_vector_at(z, i4 - 2)        = mc_blas_vector_at(z, ipn4 - i4 - 2);
					mc_blas_vector_at(z, ipn4 - i4 - 2) = temp;
					temp                                = mc_blas_vector_at(z, i4 - 1);
					mc_blas_vector_at(z, i4 - 1)        = mc_blas_vector_at(z, ipn4 - i4 - 5);
					mc_blas_vector_at(z, ipn4 - i4 - 5) = temp;
					temp                                = mc_blas_vector_at(z, i4);
					mc_blas_vector_at(z, i4)            = mc_blas_vector_at(z, ipn4 - i4 - 4);
					mc_blas_vector_at(z, ipn4 - i4 - 4) = temp;
				}
			}
		}

		dmin = -mc_fmaxf(zero, qmin - mc_sqrtf(qmin) * two * mc_sqrtf(emax));
		nbig = (n0 - i0 + 1) * 100;
		for (iwhilb = 1; iwhilb <= nbig; ++iwhilb) {
			if (i0 > n0) {
			goto F150;
		}

		mc_lapack_slasq3(i0, &n0, z, &pp, &dmin, &sigma, &desig, &qmax, &nfail, &iter, &ndiv, ieee, &ttype, &dmin1, &dmin2, &dn, &dn1, &dn2, &g, &tau);
		pp = 1 - pp;

		if (pp == 0 && n0 - i0 >= 3) {
			if (mc_blas_vector_at(z, n0 * 4) <= tol2 * qmax || mc_blas_vector_at(z, (4 * n0) - 1) <= tol2 * sigma) {
				splt   = i0 - 1;
				qmax   = mc_blas_vector_at(z, (4 * i0) - 3);
				emin   = mc_blas_vector_at(z, (4 * i0) - 1);
				oldemn = mc_blas_vector_at(z, i0 * 4);
				for (i4 = (4 * i0); i4 <= (4 * (n0 - 3)); i4 += 4) {
					if (mc_blas_vector_at(z, i4) <= tol2 * mc_blas_vector_at(z, i4 - 3) || mc_blas_vector_at(z, i4 - 1) <= tol2 * sigma) {
						mc_blas_vector_at(z, i4 - 1) = -sigma;
						splt                         = i4 / 4;
						qmax                         = zero;
						emin                         = mc_blas_vector_at(z, i4 + 3);
						oldemn                       = mc_blas_vector_at(z, i4 + 4);
					} else {
						qmax   = mc_fmaxf(qmax, mc_blas_vector_at(z, i4 + 1));
						emin   = mc_fminf(emin, mc_blas_vector_at(z, i4 - 1));
						oldemn = mc_fminf(oldemn, mc_blas_vector_at(z, i4));
					}
				}
				mc_blas_vector_at(z, (4 * n0) - 1) = emin;
				mc_blas_vector_at(z, n0 * 4)       = oldemn;
				i0                                 = splt + 1;
			}
		}
	}

	*info = 2;
	 i1   = i0;
	 n1   = n0;
	mc_unused(n1);

F145:
	tempq                              = mc_blas_vector_at(z, (4 * i0) - 3);
	mc_blas_vector_at(z, (4 * i0) - 3) = mc_blas_vector_at(z, (4 * i0) - 3) + sigma;
	for (k = i0 + 1; k <= n0; ++k) {
		tempe                             = mc_blas_vector_at(z, (4 * k) - 5);
		mc_blas_vector_at(z, (4 * k) - 5) = mc_blas_vector_at(z, (4 * k) - 5) * (tempq / mc_blas_vector_at(z, (4 * k) - 7));
		tempq                             = mc_blas_vector_at(z, (4 * k) - 3);
		mc_blas_vector_at(z, (4 * k) - 3) = mc_blas_vector_at(z, (4 * k) - 3) + sigma + tempe - mc_blas_vector_at(z, (4 * k) - 5);
	}

	if (i1 > 1) {
		n1 = i1 - 1;
		while (i1 >= 2 && mc_blas_vector_at(z, (4 * i1) - 5) >= zero) {
			i1 = i1 - 1;
		}
		if (i1 >= 1) {
			sigma = -mc_blas_vector_at(z, (4 * n1) - 1);
			goto F145;
		}
	}
	for (k = 1; k <= n; ++k) {
		mc_blas_vector_at(z, (2 * k) - 1) = mc_blas_vector_at(z, (4 * k) - 3);
		if (k < n0) {
			mc_blas_vector_at(z, k * 2) = mc_blas_vector_at(z, (4 * k) - 1);
		} else {
			mc_blas_vector_at(z, k * 2) = zero;
		}
	}
	return;

F150:;
	}

	*info = 3;
	return;

F170:
	for (k = 2; k <= n; ++k) {
		mc_blas_vector_at(z, k) = mc_blas_vector_at(z, (4 * k) - 3);
	}
	mc_lapack_slasrt('D', n, z, &iinfo);

	e = zero;
	for (k = n; k >= 1; --k) {
		e = e + mc_blas_vector_at(z, k);
	}

	mc_blas_vector_at(z, (2 * n) + 1) = trace;
	mc_blas_vector_at(z, (2 * n) + 2) = e;
	mc_blas_vector_at(z, (2 * n) + 3) = mc_cast(float, iter);
	mc_blas_vector_at(z, (2 * n) + 4) = mc_cast(float, ndiv) / mc_raise2f(mc_cast(float, n));
	mc_blas_vector_at(z, (2 * n) + 5) = nfail * hundrd / mc_cast(const float, iter);
}

#pragma mark - mc_lapack_dlasq2 -

MC_TARGET_PROC void mc_lapack_dlasq2(const int n, double * z, int * info)
{
	const double hundrd = 100.0, four = 4.0, two = 2.0, one = 1.0, half = 0.5, zero = 0.0;
	const double cbias = 1.5;

	int ieee;
	int i0, i4, iinfo, ipn4, iter, iwhila, iwhilb, k,
	    kmin, n0, nbig, ndiv, nfail, pp, splt, ttype,
	    i1, n1;
	double d, dee, deemin, desig, dmin, dmin1, dmin2, dn,
	       dn1, dn2, e, emax, emin, eps, g, oldemn, qmax,
	       qmin, s, safmin, sigma, t, tau, temp, tol,
	       tol2, trace, zmax, tempe, tempq;

	*info   = 0;
	 eps    = mc_lapack_dlamch('P');
	 safmin = mc_lapack_dlamch('S');
	 tol    = eps * hundrd;
	 tol2   = mc_raise2(tol);

	if (n < 0) {
		*info = -1;
		mc_blas_xerbla("DLASQ2", 1);
		return;
	} else if (n == 0) {
		return;
	} else if (n == 1) {
		if (mc_blas_vector_at(z, 1) < zero) {
			*info = -201;
			mc_blas_xerbla("DLASQ2", 2);
		}
		return;
	} else if (n == 2) {
		if (mc_blas_vector_at(z, 2) < zero || mc_blas_vector_at(z, 3) < zero) {
			*info = -2;
			mc_blas_xerbla("DLASQ2", 2);
			return;
		} else if (mc_blas_vector_at(z, 3) > mc_blas_vector_at(z, 1)) {
			d                       = mc_blas_vector_at(z, 3);
			mc_blas_vector_at(z, 3) = mc_blas_vector_at(z, 1);
			mc_blas_vector_at(z, 1) = d;
		}
		mc_blas_vector_at(z, 5) = mc_blas_vector_at(z, 1) + mc_blas_vector_at(z, 2) + mc_blas_vector_at(z, 3);
		if (mc_blas_vector_at(z, 2) > mc_blas_vector_at(z, 3) * tol2) {
			t = (mc_blas_vector_at(z, 1) - mc_blas_vector_at(z, 3) + mc_blas_vector_at(z, 2)) * half;
			s = mc_blas_vector_at(z, 3) * (mc_blas_vector_at(z, 2) / t);
			if (s <= t) {
				s = mc_blas_vector_at(z, 3) * (mc_blas_vector_at(z, 2) / (t * (mc_sqrt(s / t + one) + one)));
			} else {
				s = mc_blas_vector_at(z, 3) * (mc_blas_vector_at(z, 2) / (t + mc_sqrt(t) * mc_sqrt(t + s)));
			}
			t                       = mc_blas_vector_at(z, 1) + (s + mc_blas_vector_at(z, 2));
			mc_blas_vector_at(z, 3) = mc_blas_vector_at(z, 3) * (mc_blas_vector_at(z, 1) / t);
			mc_blas_vector_at(z, 1) = t;
		}
		mc_blas_vector_at(z, 2) = mc_blas_vector_at(z, 3);
		mc_blas_vector_at(z, 6) = mc_blas_vector_at(z, 2) + mc_blas_vector_at(z, 1);
		return;
	}

	mc_blas_vector_at(z, 2 * n) = zero;
	emin                        = mc_blas_vector_at(z, 2);
	qmax                        = zero;
	zmax                        = zero;
	d                           = zero;
	e                           = zero;

	for (k = 1; k <= (2 * (n - 1)); k += 2) {
		if (mc_blas_vector_at(z, k) < zero) {
			*info = -(200 + k);
			mc_blas_xerbla("DLASQ2", 2);
			return;
		} else if (mc_blas_vector_at(z, k + 1) < zero) {
			*info = -(200 + k + 1);
			mc_blas_xerbla("DLASQ2", 2);
			return;
		}
		d    = d + mc_blas_vector_at(z, k);
		e    = e + mc_blas_vector_at(z, k + 1);
		qmax = mc_fmax(qmax, mc_blas_vector_at(z, k));
		emin = mc_fmin(emin, mc_blas_vector_at(z, k + 1));
		zmax = mc_fmax(mc_maxmag(qmax, zmax), mc_blas_vector_at(z, k + 1));
	}
	if (mc_blas_vector_at(z, (2 * n) - 1) < zero) {
		*info = -(200 + ((2 * n) - 1));
		mc_blas_xerbla("DLASQ2", 2);
		return;
	}
	d    = d + mc_blas_vector_at(z, (2 * n) - 1);
	qmax = mc_fmax(qmax, mc_blas_vector_at(z, (2 * n) - 1));
	zmax = mc_fmax(qmax, zmax);
	mc_unused(zmax);

	if (e == zero) {
		for (k = 2; k <= n; ++k) {
			mc_blas_vector_at(z, k) = mc_blas_vector_at(z, (2 * k) - 1);
		}
		mc_lapack_dlasrt('D', n, z, &iinfo);
		mc_blas_vector_at(z, (2 * n) - 1) = d;
		return;
	}

	trace = d + e;

	if (trace == zero) {
		mc_blas_vector_at(z, (2 * n) - 1) = zero;
		return;
	}

	ieee = 1;

	for (k = (2 * n); k >= 2; k += -2) {
		mc_blas_vector_at(z,  2 * k     ) = zero;
		mc_blas_vector_at(z, (2 * k) - 1) = mc_blas_vector_at(z, k);
		mc_blas_vector_at(z, (2 * k) - 2) = zero;
		mc_blas_vector_at(z, (2 * k) - 3) = mc_blas_vector_at(z, k - 1);
	}

	i0 = 1;
	n0 = n;

	if (mc_blas_vector_at(z, (4 * i0) - 3) * cbias < mc_blas_vector_at(z, (4 * n0) - 3)) {
		ipn4 = 4 * (i0 + n0);
		for (i4 = (4 * i0); i4 <= (2 * (i0 + n0 - 1)); i4 += 4) {
			temp                                = mc_blas_vector_at(z, i4 - 3);
			mc_blas_vector_at(z, i4 - 3)        = mc_blas_vector_at(z, ipn4 - i4 - 3);
			mc_blas_vector_at(z, ipn4 - i4 - 3) = temp;
			temp                                = mc_blas_vector_at(z, i4 - 1);
			mc_blas_vector_at(z, i4 - 1)        = mc_blas_vector_at(z, ipn4 - i4 - 5);
			mc_blas_vector_at(z, ipn4 - i4 - 5) = temp;
		}
	}

	pp = 0;

	for (k = 1; k <= 2; ++k) {
		d = mc_blas_vector_at(z, (4 * n0) + pp - 3);
		for (i4 = ((4 * (n0 - 1)) + pp); i4 >=  ((4 * i0) + pp); i4 += -4) {
			if (mc_blas_vector_at(z, i4 - 1) <= tol2 * d) {
				mc_blas_vector_at(z, i4 - 1) = -zero;
				d                            = mc_blas_vector_at(z, i4 - 3);
			} else {
				d = mc_blas_vector_at(z, i4 - 3) * (d / (d + mc_blas_vector_at(z, i4 - 1)));
			}
		}

		emin = mc_blas_vector_at(z, (4 * i0) + pp + 1);
		d    = mc_blas_vector_at(z, (4 * i0) + pp - 3);
		for (i4 = ((4 * i0) + pp); i4 <= ((4 * (n0 - 1)) + pp); i4 += 4) {
			mc_blas_vector_at(z, i4 - (2 * pp) - 2) = d + mc_blas_vector_at(z, i4 - 1);
			if (mc_blas_vector_at(z, i4 - 1) <= tol2 * d) {
				mc_blas_vector_at(z, i4 - 1)            = -zero;
				mc_blas_vector_at(z, i4 - (2 * pp) - 2) = d;
				mc_blas_vector_at(z, i4 - (2 * pp))     = zero;
				d                                       = mc_blas_vector_at(z, i4 + 1);
			} else if (
					safmin * mc_blas_vector_at(z, i4 + 1) < mc_blas_vector_at(z, i4 - (2 * pp) - 2)
				&& safmin * mc_blas_vector_at(z, i4 - (2 * pp) - 2) < mc_blas_vector_at(z, i4 + 1)
			) {
				temp                                = mc_blas_vector_at(z, i4 + 1) / mc_blas_vector_at(z, i4 - (2 * pp) - 2);
				mc_blas_vector_at(z, i4 - (2 * pp)) = mc_blas_vector_at(z, i4 - 1) * temp;
				d                                   = d * temp;
			} else {
				mc_blas_vector_at(z, i4 - (2 * pp)) = mc_blas_vector_at(z, i4 + 1) * (mc_blas_vector_at(z, i4 - 1) / mc_blas_vector_at(z, i4 - (2 * pp) - 2));
				d                                   = mc_blas_vector_at(z, i4 + 1) * (d / mc_blas_vector_at(z, i4 - (2 * pp) - 2));
			}
			emin = mc_fmin(emin, mc_blas_vector_at(z, i4 - (2 * pp)));
		}
		mc_blas_vector_at(z, (4 * n0) - pp - 2) = d;

		qmax = mc_blas_vector_at(z, (4 * i0) - pp - 2);
		for (i4 = ((4 * i0) - pp + 2); i4 <= ((4 * n0) - pp - 2); i4 += 4) {
			qmax = mc_fmax(qmax, mc_blas_vector_at(z, i4));
		}
		pp = 1 - pp;
	}

	ttype = 0;
	dmin1 = zero;
	dmin2 = zero;
	dn    = zero;
	dn1   = zero;
	dn2   = zero;
	g     = zero;
	tau   = zero;

	iter  = 2;
	nfail = 0;
	ndiv  = 2 * (n0 - i0);

	for (iwhila = 1; iwhila <= (n + 1); ++iwhila) {
		if (n0 < 1) {
			goto F170;
		}
		desig = zero;
		if (n0 == n) {
			sigma = zero;
		} else {
			sigma = -mc_blas_vector_at(z, (4 * n0) - 1);
		}
		if (sigma < zero) {
			*info = 1;
			return;
		}
		emax = zero;
		if (n0 > i0) {
			emin = mc_fabs(mc_blas_vector_at(z, (4 * n0) - 5));
		} else {
			emin = zero;
		}
		qmin = mc_blas_vector_at(z, (4 * n0) - 3);
		qmax = qmin;
		for (i4 = (4 * n0); i4 >= 8; i4 += -4) {
			if (mc_blas_vector_at(z, i4 - 5) <= zero) {
				goto F100;
			}
			if (qmin >= emax * four) {
				qmin = mc_fmin(qmin, mc_blas_vector_at(z, i4 - 3));
				emax = mc_fmax(emax, mc_blas_vector_at(z, i4 - 5));
			}
			qmax = mc_fmax(qmax, mc_blas_vector_at(z, i4 - 7) + mc_blas_vector_at(z, i4 - 5));
			emin = mc_fmin(emin, mc_blas_vector_at(z, i4 - 5));
		}
		i4 = 4;

F100:
		i0 = i4 / 4;
		pp = 0;

		if ((n0 - i0) > 1) {
			dee    = mc_blas_vector_at(z, (4 * i0) - 3);
			deemin = dee;
			kmin   = i0;
			for (i4 = ((4 * i0) + 1); i4 <= ((4 * n0) - 3); i4 += 4) {
				dee = mc_blas_vector_at(z, i4) * (dee / (dee + mc_blas_vector_at(z, i4 - 2)));
				if (dee <= deemin) {
					deemin = dee;
					kmin = (i4 + 3) / 4;
				}
			}
			if (((kmin - i0) * 2) < n0 - kmin && deemin <= mc_blas_vector_at(z, (4 * n0) - 3) * half) {
				ipn4 = 4 * (i0 + n0);
				pp   = 2;
				for (i4 = (4 * i0); i4 <= (2 * (i0 + n0 - 1)); i4 += 4) {
					temp                                = mc_blas_vector_at(z, i4 - 3);
					mc_blas_vector_at(z, i4 - 3)        = mc_blas_vector_at(z, ipn4 - i4 - 3);
					mc_blas_vector_at(z, ipn4 - i4 - 3) = temp;
					temp                                = mc_blas_vector_at(z, i4 - 2);
					mc_blas_vector_at(z, i4 - 2)        = mc_blas_vector_at(z, ipn4 - i4 - 2);
					mc_blas_vector_at(z, ipn4 - i4 - 2) = temp;
					temp                                = mc_blas_vector_at(z, i4 - 1);
					mc_blas_vector_at(z, i4 - 1)        = mc_blas_vector_at(z, ipn4 - i4 - 5);
					mc_blas_vector_at(z, ipn4 - i4 - 5) = temp;
					temp                                = mc_blas_vector_at(z, i4);
					mc_blas_vector_at(z, i4)            = mc_blas_vector_at(z, ipn4 - i4 - 4);
					mc_blas_vector_at(z, ipn4 - i4 - 4) = temp;
				}
			}
		}

		dmin = -mc_fmax(zero, qmin - mc_sqrt(qmin) * two * mc_sqrt(emax));
		nbig = (n0 - i0 + 1) * 100;
		for (iwhilb = 1; iwhilb <= nbig; ++iwhilb) {
			if (i0 > n0) {
			goto F150;
		}

		mc_lapack_dlasq3(i0, &n0, z, &pp, &dmin, &sigma, &desig, &qmax, &nfail, &iter, &ndiv, ieee, &ttype, &dmin1, &dmin2, &dn, &dn1, &dn2, &g, &tau);
		pp = 1 - pp;

		if (pp == 0 && n0 - i0 >= 3) {
			if (mc_blas_vector_at(z, n0 * 4) <= tol2 * qmax || mc_blas_vector_at(z, (4 * n0) - 1) <= tol2 * sigma) {
				splt   = i0 - 1;
				qmax   = mc_blas_vector_at(z, (4 * i0) - 3);
				emin   = mc_blas_vector_at(z, (4 * i0) - 1);
				oldemn = mc_blas_vector_at(z, i0 * 4);
				for (i4 = (4 * i0); i4 <= (4 * (n0 - 3)); i4 += 4) {
					if (mc_blas_vector_at(z, i4) <= tol2 * mc_blas_vector_at(z, i4 - 3) || mc_blas_vector_at(z, i4 - 1) <= tol2 * sigma) {
						mc_blas_vector_at(z, i4 - 1) = -sigma;
						splt                         = i4 / 4;
						qmax                         = zero;
						emin                         = mc_blas_vector_at(z, i4 + 3);
						oldemn                       = mc_blas_vector_at(z, i4 + 4);
					} else {
						qmax   = mc_fmax(qmax, mc_blas_vector_at(z, i4 + 1));
						emin   = mc_fmin(emin, mc_blas_vector_at(z, i4 - 1));
						oldemn = mc_fmin(oldemn, mc_blas_vector_at(z, i4));
					}
				}
				mc_blas_vector_at(z, (4 * n0) - 1) = emin;
				mc_blas_vector_at(z, n0 * 4)       = oldemn;
				i0                                 = splt + 1;
			}
		}
	}

	*info = 2;
	 i1   = i0;
	 n1   = n0;
	mc_unused(n1);

F145:
	tempq                              = mc_blas_vector_at(z, (4 * i0) - 3);
	mc_blas_vector_at(z, (4 * i0) - 3) = mc_blas_vector_at(z, (4 * i0) - 3) + sigma;
	for (k = i0 + 1; k <= n0; ++k) {
		tempe                             = mc_blas_vector_at(z, (4 * k) - 5);
		mc_blas_vector_at(z, (4 * k) - 5) = mc_blas_vector_at(z, (4 * k) - 5) * (tempq / mc_blas_vector_at(z, (4 * k) - 7));
		tempq                             = mc_blas_vector_at(z, (4 * k) - 3);
		mc_blas_vector_at(z, (4 * k) - 3) = mc_blas_vector_at(z, (4 * k) - 3) + sigma + tempe - mc_blas_vector_at(z, (4 * k) - 5);
	}

	if (i1 > 1) {
		n1 = i1 - 1;
		while (i1 >= 2 && mc_blas_vector_at(z, (4 * i1) - 5) >= zero) {
			i1 = i1 - 1;
		}
		if (i1 >= 1) {
			sigma = -mc_blas_vector_at(z, (4 * n1) - 1);
			goto F145;
		}
	}
	for (k = 1; k <= n; ++k) {
		mc_blas_vector_at(z, (2 * k) - 1) = mc_blas_vector_at(z, (4 * k) - 3);
		if (k < n0) {
			mc_blas_vector_at(z, k * 2) = mc_blas_vector_at(z, (4 * k) - 1);
		} else {
			mc_blas_vector_at(z, k * 2) = zero;
		}
	}
	return;

F150:;
	}

	*info = 3;
	return;

F170:
	for (k = 2; k <= n; ++k) {
		mc_blas_vector_at(z, k) = mc_blas_vector_at(z, (4 * k) - 3);
	}
	mc_lapack_dlasrt('D', n, z, &iinfo);

	e = zero;
	for (k = n; k >= 1; --k) {
		e = e + mc_blas_vector_at(z, k);
	}

	mc_blas_vector_at(z, (2 * n) + 1) = trace;
	mc_blas_vector_at(z, (2 * n) + 2) = e;
	mc_blas_vector_at(z, (2 * n) + 3) = mc_cast(double, iter);
	mc_blas_vector_at(z, (2 * n) + 4) = mc_cast(double, ndiv) / mc_raise2(mc_cast(double, n));
	mc_blas_vector_at(z, (2 * n) + 5) = nfail * hundrd / mc_cast(const double, iter);
}

#pragma mark - mc_lapack_llasq2 -

MC_TARGET_PROC void mc_lapack_llasq2(const int n, long double * z, int * info)
{
	const long double hundrd = 100.0L, four = 4.0L, two = 2.0L, one = 1.0L, half = 0.5L, zero = 0.0L;
	const long double cbias = 1.5L;

	int ieee;
	int i0, i4, iinfo, ipn4, iter, iwhila, iwhilb, k,
	    kmin, n0, nbig, ndiv, nfail, pp, splt, ttype,
	    i1, n1;
	long double d, dee, deemin, desig, dmin, dmin1, dmin2, dn,
	            dn1, dn2, e, emax, emin, eps, g, oldemn, qmax,
	            qmin, s, safmin, sigma, t, tau, temp, tol,
	            tol2, trace, zmax, tempe, tempq;

	*info   = 0;
	 eps    = mc_lapack_llamch('P');
	 safmin = mc_lapack_llamch('S');
	 tol    = eps * hundrd;
	 tol2   = mc_raise2l(tol);

	if (n < 0) {
		*info = -1;
		mc_blas_xerbla("LLASQ2", 1);
		return;
	} else if (n == 0) {
		return;
	} else if (n == 1) {
		if (mc_blas_vector_at(z, 1) < zero) {
			*info = -201;
			mc_blas_xerbla("LLASQ2", 2);
		}
		return;
	} else if (n == 2) {
		if (mc_blas_vector_at(z, 2) < zero || mc_blas_vector_at(z, 3) < zero) {
			*info = -2;
			mc_blas_xerbla("LLASQ2", 2);
			return;
		} else if (mc_blas_vector_at(z, 3) > mc_blas_vector_at(z, 1)) {
			d                       = mc_blas_vector_at(z, 3);
			mc_blas_vector_at(z, 3) = mc_blas_vector_at(z, 1);
			mc_blas_vector_at(z, 1) = d;
		}
		mc_blas_vector_at(z, 5) = mc_blas_vector_at(z, 1) + mc_blas_vector_at(z, 2) + mc_blas_vector_at(z, 3);
		if (mc_blas_vector_at(z, 2) > mc_blas_vector_at(z, 3) * tol2) {
			t = (mc_blas_vector_at(z, 1) - mc_blas_vector_at(z, 3) + mc_blas_vector_at(z, 2)) * half;
			s = mc_blas_vector_at(z, 3) * (mc_blas_vector_at(z, 2) / t);
			if (s <= t) {
				s = mc_blas_vector_at(z, 3) * (mc_blas_vector_at(z, 2) / (t * (mc_sqrtl(s / t + one) + one)));
			} else {
				s = mc_blas_vector_at(z, 3) * (mc_blas_vector_at(z, 2) / (t + mc_sqrtl(t) * mc_sqrtl(t + s)));
			}
			t                       = mc_blas_vector_at(z, 1) + (s + mc_blas_vector_at(z, 2));
			mc_blas_vector_at(z, 3) = mc_blas_vector_at(z, 3) * (mc_blas_vector_at(z, 1) / t);
			mc_blas_vector_at(z, 1) = t;
		}
		mc_blas_vector_at(z, 2) = mc_blas_vector_at(z, 3);
		mc_blas_vector_at(z, 6) = mc_blas_vector_at(z, 2) + mc_blas_vector_at(z, 1);
		return;
	}

	mc_blas_vector_at(z, 2 * n) = zero;
	emin                        = mc_blas_vector_at(z, 2);
	qmax                        = zero;
	zmax                        = zero;
	d                           = zero;
	e                           = zero;

	for (k = 1; k <= (2 * (n - 1)); k += 2) {
		if (mc_blas_vector_at(z, k) < zero) {
			*info = -(200 + k);
			mc_blas_xerbla("LLASQ2", 2);
			return;
		} else if (mc_blas_vector_at(z, k + 1) < zero) {
			*info = -(200 + k + 1);
			mc_blas_xerbla("LLASQ2", 2);
			return;
		}
		d    = d + mc_blas_vector_at(z, k);
		e    = e + mc_blas_vector_at(z, k + 1);
		qmax = mc_fmaxl(qmax, mc_blas_vector_at(z, k));
		emin = mc_fminl(emin, mc_blas_vector_at(z, k + 1));
		zmax = mc_fmaxl(mc_maxmag(qmax, zmax), mc_blas_vector_at(z, k + 1));
	}
	if (mc_blas_vector_at(z, (2 * n) - 1) < zero) {
		*info = -(200 + ((2 * n) - 1));
		mc_blas_xerbla("LLASQ2", 2);
		return;
	}
	d    = d + mc_blas_vector_at(z, (2 * n) - 1);
	qmax = mc_fmaxl(qmax, mc_blas_vector_at(z, (2 * n) - 1));
	zmax = mc_fmaxl(qmax, zmax);
	mc_unused(zmax);

	if (e == zero) {
		for (k = 2; k <= n; ++k) {
			mc_blas_vector_at(z, k) = mc_blas_vector_at(z, (2 * k) - 1);
		}
		mc_lapack_llasrt('D', n, z, &iinfo);
		mc_blas_vector_at(z, (2 * n) - 1) = d;
		return;
	}

	trace = d + e;

	if (trace == zero) {
		mc_blas_vector_at(z, (2 * n) - 1) = zero;
		return;
	}

	ieee = 1;

	for (k = (2 * n); k >= 2; k += -2) {
		mc_blas_vector_at(z,  2 * k     ) = zero;
		mc_blas_vector_at(z, (2 * k) - 1) = mc_blas_vector_at(z, k);
		mc_blas_vector_at(z, (2 * k) - 2) = zero;
		mc_blas_vector_at(z, (2 * k) - 3) = mc_blas_vector_at(z, k - 1);
	}

	i0 = 1;
	n0 = n;

	if (mc_blas_vector_at(z, (4 * i0) - 3) * cbias < mc_blas_vector_at(z, (4 * n0) - 3)) {
		ipn4 = 4 * (i0 + n0);
		for (i4 = (4 * i0); i4 <= (2 * (i0 + n0 - 1)); i4 += 4) {
			temp                                = mc_blas_vector_at(z, i4 - 3);
			mc_blas_vector_at(z, i4 - 3)        = mc_blas_vector_at(z, ipn4 - i4 - 3);
			mc_blas_vector_at(z, ipn4 - i4 - 3) = temp;
			temp                                = mc_blas_vector_at(z, i4 - 1);
			mc_blas_vector_at(z, i4 - 1)        = mc_blas_vector_at(z, ipn4 - i4 - 5);
			mc_blas_vector_at(z, ipn4 - i4 - 5) = temp;
		}
	}

	pp = 0;

	for (k = 1; k <= 2; ++k) {
		d = mc_blas_vector_at(z, (4 * n0) + pp - 3);
		for (i4 = ((4 * (n0 - 1)) + pp); i4 >=  ((4 * i0) + pp); i4 += -4) {
			if (mc_blas_vector_at(z, i4 - 1) <= tol2 * d) {
				mc_blas_vector_at(z, i4 - 1) = -zero;
				d                            = mc_blas_vector_at(z, i4 - 3);
			} else {
				d = mc_blas_vector_at(z, i4 - 3) * (d / (d + mc_blas_vector_at(z, i4 - 1)));
			}
		}

		emin = mc_blas_vector_at(z, (4 * i0) + pp + 1);
		d    = mc_blas_vector_at(z, (4 * i0) + pp - 3);
		for (i4 = ((4 * i0) + pp); i4 <= ((4 * (n0 - 1)) + pp); i4 += 4) {
			mc_blas_vector_at(z, i4 - (2 * pp) - 2) = d + mc_blas_vector_at(z, i4 - 1);
			if (mc_blas_vector_at(z, i4 - 1) <= tol2 * d) {
				mc_blas_vector_at(z, i4 - 1)            = -zero;
				mc_blas_vector_at(z, i4 - (2 * pp) - 2) = d;
				mc_blas_vector_at(z, i4 - (2 * pp))     = zero;
				d                                       = mc_blas_vector_at(z, i4 + 1);
			} else if (
					safmin * mc_blas_vector_at(z, i4 + 1) < mc_blas_vector_at(z, i4 - (2 * pp) - 2)
				&& safmin * mc_blas_vector_at(z, i4 - (2 * pp) - 2) < mc_blas_vector_at(z, i4 + 1)
			) {
				temp                                = mc_blas_vector_at(z, i4 + 1) / mc_blas_vector_at(z, i4 - (2 * pp) - 2);
				mc_blas_vector_at(z, i4 - (2 * pp)) = mc_blas_vector_at(z, i4 - 1) * temp;
				d                                   = d * temp;
			} else {
				mc_blas_vector_at(z, i4 - (2 * pp)) = mc_blas_vector_at(z, i4 + 1) * (mc_blas_vector_at(z, i4 - 1) / mc_blas_vector_at(z, i4 - (2 * pp) - 2));
				d                                   = mc_blas_vector_at(z, i4 + 1) * (d / mc_blas_vector_at(z, i4 - (2 * pp) - 2));
			}
			emin = mc_fminl(emin, mc_blas_vector_at(z, i4 - (2 * pp)));
		}
		mc_blas_vector_at(z, (4 * n0) - pp - 2) = d;

		qmax = mc_blas_vector_at(z, (4 * i0) - pp - 2);
		for (i4 = ((4 * i0) - pp + 2); i4 <= ((4 * n0) - pp - 2); i4 += 4) {
			qmax = mc_fmaxl(qmax, mc_blas_vector_at(z, i4));
		}
		pp = 1 - pp;
	}

	ttype = 0;
	dmin1 = zero;
	dmin2 = zero;
	dn    = zero;
	dn1   = zero;
	dn2   = zero;
	g     = zero;
	tau   = zero;

	iter  = 2;
	nfail = 0;
	ndiv  = 2 * (n0 - i0);

	for (iwhila = 1; iwhila <= (n + 1); ++iwhila) {
		if (n0 < 1) {
			goto F170;
		}
		desig = zero;
		if (n0 == n) {
			sigma = zero;
		} else {
			sigma = -mc_blas_vector_at(z, (4 * n0) - 1);
		}
		if (sigma < zero) {
			*info = 1;
			return;
		}
		emax = zero;
		if (n0 > i0) {
			emin = mc_fabsl(mc_blas_vector_at(z, (4 * n0) - 5));
		} else {
			emin = zero;
		}
		qmin = mc_blas_vector_at(z, (4 * n0) - 3);
		qmax = qmin;
		for (i4 = (4 * n0); i4 >= 8; i4 += -4) {
			if (mc_blas_vector_at(z, i4 - 5) <= zero) {
				goto F100;
			}
			if (qmin >= emax * four) {
				qmin = mc_fminl(qmin, mc_blas_vector_at(z, i4 - 3));
				emax = mc_fmaxl(emax, mc_blas_vector_at(z, i4 - 5));
			}
			qmax = mc_fmaxl(qmax, mc_blas_vector_at(z, i4 - 7) + mc_blas_vector_at(z, i4 - 5));
			emin = mc_fminl(emin, mc_blas_vector_at(z, i4 - 5));
		}
		i4 = 4;

F100:
		i0 = i4 / 4;
		pp = 0;

		if ((n0 - i0) > 1) {
			dee    = mc_blas_vector_at(z, (4 * i0) - 3);
			deemin = dee;
			kmin   = i0;
			for (i4 = ((4 * i0) + 1); i4 <= ((4 * n0) - 3); i4 += 4) {
				dee = mc_blas_vector_at(z, i4) * (dee / (dee + mc_blas_vector_at(z, i4 - 2)));
				if (dee <= deemin) {
					deemin = dee;
					kmin = (i4 + 3) / 4;
				}
			}
			if (((kmin - i0) * 2) < n0 - kmin && deemin <= mc_blas_vector_at(z, (4 * n0) - 3) * half) {
				ipn4 = 4 * (i0 + n0);
				pp   = 2;
				for (i4 = (4 * i0); i4 <= (2 * (i0 + n0 - 1)); i4 += 4) {
					temp                                = mc_blas_vector_at(z, i4 - 3);
					mc_blas_vector_at(z, i4 - 3)        = mc_blas_vector_at(z, ipn4 - i4 - 3);
					mc_blas_vector_at(z, ipn4 - i4 - 3) = temp;
					temp                                = mc_blas_vector_at(z, i4 - 2);
					mc_blas_vector_at(z, i4 - 2)        = mc_blas_vector_at(z, ipn4 - i4 - 2);
					mc_blas_vector_at(z, ipn4 - i4 - 2) = temp;
					temp                                = mc_blas_vector_at(z, i4 - 1);
					mc_blas_vector_at(z, i4 - 1)        = mc_blas_vector_at(z, ipn4 - i4 - 5);
					mc_blas_vector_at(z, ipn4 - i4 - 5) = temp;
					temp                                = mc_blas_vector_at(z, i4);
					mc_blas_vector_at(z, i4)            = mc_blas_vector_at(z, ipn4 - i4 - 4);
					mc_blas_vector_at(z, ipn4 - i4 - 4) = temp;
				}
			}
		}

		dmin = -mc_fmaxl(zero, qmin - mc_sqrtl(qmin) * two * mc_sqrtl(emax));
		nbig = (n0 - i0 + 1) * 100;
		for (iwhilb = 1; iwhilb <= nbig; ++iwhilb) {
			if (i0 > n0) {
			goto F150;
		}

		mc_lapack_llasq3(i0, &n0, z, &pp, &dmin, &sigma, &desig, &qmax, &nfail, &iter, &ndiv, ieee, &ttype, &dmin1, &dmin2, &dn, &dn1, &dn2, &g, &tau);
		pp = 1 - pp;

		if (pp == 0 && n0 - i0 >= 3) {
			if (mc_blas_vector_at(z, n0 * 4) <= tol2 * qmax || mc_blas_vector_at(z, (4 * n0) - 1) <= tol2 * sigma) {
				splt   = i0 - 1;
				qmax   = mc_blas_vector_at(z, (4 * i0) - 3);
				emin   = mc_blas_vector_at(z, (4 * i0) - 1);
				oldemn = mc_blas_vector_at(z, i0 * 4);
				for (i4 = (4 * i0); i4 <= (4 * (n0 - 3)); i4 += 4) {
					if (mc_blas_vector_at(z, i4) <= tol2 * mc_blas_vector_at(z, i4 - 3) || mc_blas_vector_at(z, i4 - 1) <= tol2 * sigma) {
						mc_blas_vector_at(z, i4 - 1) = -sigma;
						splt                         = i4 / 4;
						qmax                         = zero;
						emin                         = mc_blas_vector_at(z, i4 + 3);
						oldemn                       = mc_blas_vector_at(z, i4 + 4);
					} else {
						qmax   = mc_fmaxl(qmax, mc_blas_vector_at(z, i4 + 1));
						emin   = mc_fminl(emin, mc_blas_vector_at(z, i4 - 1));
						oldemn = mc_fminl(oldemn, mc_blas_vector_at(z, i4));
					}
				}
				mc_blas_vector_at(z, (4 * n0) - 1) = emin;
				mc_blas_vector_at(z, n0 * 4)       = oldemn;
				i0                                 = splt + 1;
			}
		}
	}

	*info = 2;
	 i1   = i0;
	 n1   = n0;
	mc_unused(n1);

F145:
	tempq                              = mc_blas_vector_at(z, (4 * i0) - 3);
	mc_blas_vector_at(z, (4 * i0) - 3) = mc_blas_vector_at(z, (4 * i0) - 3) + sigma;
	for (k = i0 + 1; k <= n0; ++k) {
		tempe                             = mc_blas_vector_at(z, (4 * k) - 5);
		mc_blas_vector_at(z, (4 * k) - 5) = mc_blas_vector_at(z, (4 * k) - 5) * (tempq / mc_blas_vector_at(z, (4 * k) - 7));
		tempq                             = mc_blas_vector_at(z, (4 * k) - 3);
		mc_blas_vector_at(z, (4 * k) - 3) = mc_blas_vector_at(z, (4 * k) - 3) + sigma + tempe - mc_blas_vector_at(z, (4 * k) - 5);
	}

	if (i1 > 1) {
		n1 = i1 - 1;
		while (i1 >= 2 && mc_blas_vector_at(z, (4 * i1) - 5) >= zero) {
			i1 = i1 - 1;
		}
		if (i1 >= 1) {
			sigma = -mc_blas_vector_at(z, (4 * n1) - 1);
			goto F145;
		}
	}
	for (k = 1; k <= n; ++k) {
		mc_blas_vector_at(z, (2 * k) - 1) = mc_blas_vector_at(z, (4 * k) - 3);
		if (k < n0) {
			mc_blas_vector_at(z, k * 2) = mc_blas_vector_at(z, (4 * k) - 1);
		} else {
			mc_blas_vector_at(z, k * 2) = zero;
		}
	}
	return;

F150:;
	}

	*info = 3;
	return;

F170:
	for (k = 2; k <= n; ++k) {
		mc_blas_vector_at(z, k) = mc_blas_vector_at(z, (4 * k) - 3);
	}
	mc_lapack_llasrt('D', n, z, &iinfo);

	e = zero;
	for (k = n; k >= 1; --k) {
		e = e + mc_blas_vector_at(z, k);
	}

	mc_blas_vector_at(z, (2 * n) + 1) = trace;
	mc_blas_vector_at(z, (2 * n) + 2) = e;
	mc_blas_vector_at(z, (2 * n) + 3) = mc_cast(long double, iter);
	mc_blas_vector_at(z, (2 * n) + 4) = mc_cast(long double, ndiv) / mc_raise2l(mc_cast(long double, n));
	mc_blas_vector_at(z, (2 * n) + 5) = nfail * hundrd / mc_cast(const long double, iter);
}

#endif /* !MC_LAPACKE_LASQ2_H */

/* EOF */
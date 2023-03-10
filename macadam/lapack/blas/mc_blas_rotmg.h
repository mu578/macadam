//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_rotmg.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?rotmg constructs a Gentleman's modified Given's plane rotation.
 *
 * \synopsis
 *    real-floating ?rotmg(d1, d2, x1, y1, param)
 *    real-floating d1, d2, x1, y1
 *    real-floating param(5)
 *
 * \purpose
 *    ?rotmg computes the parameters for a modified Givens rotation. Construct Gentleman's
 *    modified a Given's plane rotation that will annihilate an element of a vector:
 *       flag=-1:
 *       H[2x2] = | h11 | h12 |
 *                | h21 | h22 |
 *       flag=0:
 *       H[2x2] = |  1  | h12 |
 *                | h21 |  1  |
 *       flag=1:
 *       H[2x2] = | h11 |  1  |
 *                | -1  | h22 |
 *       flag=-2:
 *       H[2x2] = |  1  |  0  |
 *                |  0  |  1  |
 *
 * \parameters
 *    [out] d1    - real-floating. The first diagonal entry in the H matrix. Overwritten to
 *    reflect the effect of the transformation.
 *
 *    [out] d2    - real-floating. The second diagonal entry in the H matrix. Overwritten to
 *    reflect the effect of the transformation.
 *
 *    [out] x1    - real-floating. The first element of the vector to which the H matrix is
 *    applied.  to reflect the effect of the transformation.
 *
 *    [in] y1     - real-floating. The second element of the vector to which the H matrix is
 *    applied.
 *
 *    [out] param - real-floating array of size 5.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 */

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>

#ifndef MC_BLAS_ROTMG_H
#define MC_BLAS_ROTMG_H

#pragma mark - mc_blas_srotmg -

MC_TARGET_FUNC void mc_blas_srotmg(float * d1, float * d2, float * x1, const float y1, float param[5])
{
	const float zero = 0.0f, one = 1.0f, two = 2.0f;

	const float gam    = +4096.00000000000000000000000000000000000E+00f;
	const float gamsq  = +16777200.0000000000000000000000000000000E+00f;
	const float rgamsq = 1.0f / gamsq;

	float flag, h11, h12, h21, h22, p1, p2, q1, q2, temp, u;

	if ((*d1) < zero) {
		flag  = -one;
		h11   = zero;
		h12   = zero;
		h21   = zero;
		h22   = zero;

		(*d1) = zero;
		(*d2) = zero;
		(*x1) = zero;
	} else {
		p2 = *d2 * y1;
		if (p2 == zero) {
			flag                        = -two;
			mc_blas_vector_at(param, 1) = flag;
			return;
		}
		p1 = (*d1) * (*x1);
		q2 = p2 * y1;
		q1 = p1 * (*x1);
		if (mc_fabsf(q1) > mc_fabsf(q2)) {
			h21 = -y1 / (*x1);
			h12 = p2 / p1;
			u   = one - h12 * h21;
			if (u > zero) {
				flag  = zero;
				(*d1) = (*d1) / u;
				(*d2) = (*d2) / u;
				(*x1) = (*x1) * u;
			}
		} else {
			if (q2 < zero) {
				flag  = -one;
				h11   = zero;
				h12   = zero;
				h21   = zero;
				h22   = zero;

				(*d1) = zero;
				(*d2) = zero;
				(*x1) = zero;
			} else {
				flag  = one;
				h11   = p1 / p2;
				h22   = (*x1) / y1;
				u     = one + h11 * h22;
				temp  = (*d2) / u;
				(*d2) = (*d1) / u;
				(*d1) = temp;
				(*x1) = y1 * u;
			}
		}
		if ((*d1) != zero) {
			while ((*d1) <= rgamsq || (*d1) >= gamsq) {
				if (flag == zero) {
					h11  = one;
					h22  = one;
					flag = -one;
				} else {
					h21  = -one;
					h12  = one;
					flag = -one;
				}
				if ((*d1) <= rgamsq) {
					(*d1) = (*d1) * mc_raise2f(gam);
					(*x1) = (*x1) / gam;
					h11   = h11 / gam;
					h12   = h12 / gam;
				} else {
					(*d1) = (*d1) / mc_raise2f(gam);
					(*x1) = (*x1) * gam;
					h11   = h11 * gam;
					h12   = h12 * gam;
				}
			}
		}
		if ((*d2) != zero) {
			while (mc_fabsf(*d2) <= rgamsq || mc_fabsf(*d2) >= gamsq) {
				if (flag == zero) {
					h11  = one;
					h22  = one;
					flag = -one;
				} else {
					h21  = -one;
					h12  = one;
					flag = -one;
				}
				if (mc_fabsf(*d2) <= rgamsq) {
					(*d2) = (*d2) * mc_raise2f(gam);
					h21   = h21 / gam;
					h22   = h22 / gam;
				} else {
					(*d2) = (*d2) / mc_raise2f(gam);
					h21   = h21 * gam;
					h22   = h22 * gam;
				}
			}
		}
	}
	if (flag < zero) {
		mc_blas_vector_at(param, 2) = h11;
		mc_blas_vector_at(param, 3) = h21;
		mc_blas_vector_at(param, 4) = h12;
		mc_blas_vector_at(param, 5) = h22;
	} else if (flag == zero) {
		mc_blas_vector_at(param, 3) = h21;
		mc_blas_vector_at(param, 4) = h12;
	} else {
		mc_blas_vector_at(param, 2) = h11;
		mc_blas_vector_at(param, 5) = h22;
	}
	mc_blas_vector_at(param, 1) = flag;

#	if MC_TARGET_BLAS_USE_CLAYOUT
	mcswap_var(temp, mc_blas_vector_at(param, 3), mc_blas_vector_at(param, 4));
#	endif
}

#pragma mark - mc_blas_drotmg -

MC_TARGET_FUNC void mc_blas_drotmg(double * d1, double * d2, double * x1, const double y1, double param[5])
{
	const double zero = 0.0, one = 1.0, two = 2.0;

	const double gam    = +4096.0000000000000000000000000000000000000E+00;
	const double gamsq  = +16777200.000000000000000000000000000000000E+00;
	const double rgamsq = 1.0 / gamsq;

	double flag, h11, h12, h21, h22, p1, p2, q1, q2, temp, u;

	if ((*d1) < zero) {
		flag  = -one;
		h11   = zero;
		h12   = zero;
		h21   = zero;
		h22   = zero;

		(*d1) = zero;
		(*d2) = zero;
		(*x1) = zero;
	} else {
		p2 = *d2 * y1;
		if (p2 == zero) {
			flag                        = -two;
			mc_blas_vector_at(param, 1) = flag;
			return;
		}
		p1 = (*d1) * (*x1);
		q2 = p2 * y1;
		q1 = p1 * (*x1);
		if (mc_fabs(q1) > mc_fabs(q2)) {
			h21 = -y1 / (*x1);
			h12 = p2 / p1;
			u   = one - h12 * h21;
			if (u > zero) {
				flag  = zero;
				(*d1) = (*d1) / u;
				(*d2) = (*d2) / u;
				(*x1) = (*x1) * u;
			}
		} else {
			if (q2 < zero) {
				flag  = -one;
				h11   = zero;
				h12   = zero;
				h21   = zero;
				h22   = zero;

				(*d1) = zero;
				(*d2) = zero;
				(*x1) = zero;
			} else {
				flag  = one;
				h11   = p1 / p2;
				h22   = (*x1) / y1;
				u     = one + h11 * h22;
				temp  = (*d2) / u;
				(*d2) = (*d1) / u;
				(*d1) = temp;
				(*x1) = y1 * u;
			}
		}
		if ((*d1) != zero) {
			while ((*d1) <= rgamsq || (*d1) >= gamsq) {
				if (flag == zero) {
					h11  = one;
					h22  = one;
					flag = -one;
				} else {
					h21  = -one;
					h12  = one;
					flag = -one;
				}
				if ((*d1) <= rgamsq) {
					(*d1) = (*d1) * mc_raise2(gam);
					(*x1) = (*x1) / gam;
					h11   = h11 / gam;
					h12   = h12 / gam;
				} else {
					(*d1) = (*d1) / mc_raise2(gam);
					(*x1) = (*x1) * gam;
					h11   = h11 * gam;
					h12   = h12 * gam;
				}
			}
		}
		if ((*d2) != zero) {
			while (mc_fabs(*d2) <= rgamsq || mc_fabs(*d2) >= gamsq) {
				if (flag == zero) {
					h11  = one;
					h22  = one;
					flag = -one;
				} else {
					h21  = -one;
					h12  = one;
					flag = -one;
				}
				if (mc_fabs(*d2) <= rgamsq) {
					(*d2) = (*d2) * mc_raise2(gam);
					h21   = h21 / gam;
					h22   = h22 / gam;
				} else {
					(*d2) = (*d2) / mc_raise2(gam);
					h21   = h21 * gam;
					h22   = h22 * gam;
				}
			}
		}
	}
	if (flag < zero) {
		mc_blas_vector_at(param, 2) = h11;
		mc_blas_vector_at(param, 3) = h21;
		mc_blas_vector_at(param, 4) = h12;
		mc_blas_vector_at(param, 5) = h22;
	} else if (flag == zero) {
		mc_blas_vector_at(param, 3) = h21;
		mc_blas_vector_at(param, 4) = h12;
	} else {
		mc_blas_vector_at(param, 2) = h11;
		mc_blas_vector_at(param, 5) = h22;
	}
	mc_blas_vector_at(param, 1) = flag;

#	if MC_TARGET_BLAS_USE_CLAYOUT
	mcswap_var(temp, mc_blas_vector_at(param, 3), mc_blas_vector_at(param, 4));
#	endif
}

#pragma mark - mc_blas_lrotmg -

MC_TARGET_FUNC void mc_blas_lrotmg(long double * d1, long double * d2, long double * x1, const long double y1, long double param[5])
{
	const long double zero = 0.0L, one = 1.0L, two = 2.0L;

#	if MC_TARGET_HAVE_LONG_DOUBLE
	const long double gam    = +4096.000000000000000000000000000000000000000000000000000000000000E+00L;
	const long double gamsq  = +16777200.00000000000000000000000000000000000000000000000000000000E+00L;
	const long double rgamsq = 1.0L / gamsq;
#	else
	const long double gam    = +4096.0000000000000000000000000000000000000E+00L;
	const long double gamsq  = +16777200.000000000000000000000000000000000E+00L;
	const long double rgamsq = 1.0L / gamsq;
#	endif

	long double flag, h11, h12, h21, h22, p1, p2, q1, q2, temp, u;

	if ((*d1) < zero) {
		flag  = -one;
		h11   = zero;
		h12   = zero;
		h21   = zero;
		h22   = zero;

		(*d1) = zero;
		(*d2) = zero;
		(*x1) = zero;
	} else {
		p2 = *d2 * y1;
		if (p2 == zero) {
			flag                        = -two;
			mc_blas_vector_at(param, 1) = flag;
			return;
		}
		p1 = (*d1) * (*x1);
		q2 = p2 * y1;
		q1 = p1 * (*x1);
		if (mc_fabsl(q1) > mc_fabsl(q2)) {
			h21 = -y1 / (*x1);
			h12 = p2 / p1;
			u   = one - h12 * h21;
			if (u > zero) {
				flag  = zero;
				(*d1) = (*d1) / u;
				(*d2) = (*d2) / u;
				(*x1) = (*x1) * u;
			}
		} else {
			if (q2 < zero) {
				flag  = -one;
				h11   = zero;
				h12   = zero;
				h21   = zero;
				h22   = zero;

				(*d1) = zero;
				(*d2) = zero;
				(*x1) = zero;
			} else {
				flag  = one;
				h11   = p1 / p2;
				h22   = (*x1) / y1;
				u     = one + h11 * h22;
				temp  = (*d2) / u;
				(*d2) = (*d1) / u;
				(*d1) = temp;
				(*x1) = y1 * u;
			}
		}
		if ((*d1) != zero) {
			while ((*d1) <= rgamsq || (*d1) >= gamsq) {
				if (flag == zero) {
					h11  = one;
					h22  = one;
					flag = -one;
				} else {
					h21  = -one;
					h12  = one;
					flag = -one;
				}
				if ((*d1) <= rgamsq) {
					(*d1) = (*d1) * mc_raise2l(gam);
					(*x1) = (*x1) / gam;
					h11   = h11 / gam;
					h12   = h12 / gam;
				} else {
					(*d1) = (*d1) / mc_raise2l(gam);
					(*x1) = (*x1) * gam;
					h11   = h11 * gam;
					h12   = h12 * gam;
				}
			}
		}
		if ((*d2) != zero) {
			while (mc_fabsl(*d2) <= rgamsq || mc_fabsl(*d2) >= gamsq) {
				if (flag == zero) {
					h11  = one;
					h22  = one;
					flag = -one;
				} else {
					h21  = -one;
					h12  = one;
					flag = -one;
				}
				if (mc_fabsl(*d2) <= rgamsq) {
					(*d2) = (*d2) * mc_raise2l(gam);
					h21   = h21 / gam;
					h22   = h22 / gam;
				} else {
					(*d2) = (*d2) / mc_raise2l(gam);
					h21   = h21 * gam;
					h22   = h22 * gam;
				}
			}
		}
	}
	if (flag < zero) {
		mc_blas_vector_at(param, 2) = h11;
		mc_blas_vector_at(param, 3) = h21;
		mc_blas_vector_at(param, 4) = h12;
		mc_blas_vector_at(param, 5) = h22;
	} else if (flag == zero) {
		mc_blas_vector_at(param, 3) = h21;
		mc_blas_vector_at(param, 4) = h12;
	} else {
		mc_blas_vector_at(param, 2) = h11;
		mc_blas_vector_at(param, 5) = h22;
	}
	mc_blas_vector_at(param, 1) = flag;

#	if MC_TARGET_BLAS_USE_CLAYOUT
	mcswap_var(temp, mc_blas_vector_at(param, 3), mc_blas_vector_at(param, 4));
#	endif
}

#endif /* !MC_BLAS_ROTMG_H */

/* EOF */
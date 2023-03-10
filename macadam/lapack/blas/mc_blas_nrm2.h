//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_nrm2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    ?nrm2 returns the euclidean norm of a vector.
 *
 * \synopsis
 *    real-floating ?nrm2(n, x, incx)
 *    int           incx, n
 *    real-floating x(*)
 *
 * \purpose
 *    ?nrm2 returns the euclidean norm of a vector: norm2=sqrt(x'*x).
 *
 * \parameters
 *    [in] n     - int. Specifies the number of elements in the input vector `x`.
 *    [in] x     - real-floating array of size at least (1+(n-1)*abs(incx)).
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Sven Hammarling, Nag Ltd.
 */

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>

#ifndef MC_BLAS_NRM2_H
#define MC_BLAS_NRM2_H

#pragma mark - mc_blas_snrm2 -

MC_TARGET_FUNC float mc_blas_snrm2(const int n, const float * x, const int incx)
{
	const float one = 1.0f, zero = 0.0f;

	float absxi, norm, scale, ssq;
	int ix;

	if (n < 1 || incx < 1) {
		norm = zero;
	} else if (n == 1) {
		norm = mc_fabsf(mc_blas_vector_at(x, 1));
	} else {
		scale = zero;
		ssq   = one;
		for (ix = 1; incx < 0 ? ix >= (n - 1) * incx + 1 : ix <= (n - 1) * incx + 1; ix += incx) {
			absxi = mc_blas_vector_at(x, ix);
			if (absxi != zero) {
				absxi = mc_fabsf(absxi);
				if (scale < absxi) {
					ssq   = one + (ssq * mc_raise2f(scale / absxi));
					scale = absxi;
				} else {
					ssq = ssq + mc_raise2f(absxi / scale);
				}
			}
		}
		norm = scale * mc_sqrtf(ssq);
	}
	return norm;
}

#pragma mark - mc_blas_dsnrm2 -

MC_TARGET_FUNC double mc_blas_dsnrm2(const int n, const float * x, const int incx)
{
	const double one = 1.0, zero = 0.0;

	double absxi, norm, scale, ssq;
	int ix;

	if (n < 1 || incx < 1) {
		norm = zero;
	} else if (n == 1) {
		norm = mc_fabs(mc_cast(double, mc_blas_vector_at(x, 1)));
	} else {
		scale = zero;
		ssq   = one;
		for (ix = 1; incx < 0 ? ix >= (n - 1) * incx + 1 : ix <= (n - 1) * incx + 1; ix += incx) {
			absxi = mc_cast(double, mc_blas_vector_at(x, ix));
			if (absxi != zero) {
				absxi = mc_fabs(absxi);
				if (scale < absxi) {
					ssq   = one + (ssq * mc_raise2(scale / absxi));
					scale = absxi;
				} else {
					ssq = ssq + mc_raise2(absxi / scale);
				}
			}
		}
		norm = scale * mc_sqrt(ssq);
	}
	return norm;
}

#pragma mark - mc_blas_dnrm2 -

MC_TARGET_FUNC double mc_blas_dnrm2(const int n, const double * x, const int incx)
{
	const double one = 1.0, zero = 0.0;

	double absxi, norm, scale, ssq;
	int ix;

	if (n < 1 || incx < 1) {
		norm = zero;
	} else if (n == 1) {
		norm = mc_fabs(mc_blas_vector_at(x, 1));
	} else {
		scale = zero;
		ssq   = one;
		for (ix = 1; incx < 0 ? ix >= (n - 1) * incx + 1 : ix <= (n - 1) * incx + 1; ix += incx) {
			absxi = mc_blas_vector_at(x, ix);
			if (absxi != zero) {
				absxi = mc_fabs(absxi);
				if (scale < absxi) {
					ssq   = one + (ssq * mc_raise2(scale / absxi));
					scale = absxi;
				} else {
					ssq = ssq + mc_raise2(absxi / scale);
				}
			}
		}
		norm = scale * mc_sqrt(ssq);
	}
	return norm;
}

#pragma mark - mc_blas_lnrm2 -

MC_TARGET_FUNC long double mc_blas_lnrm2(const int n, const long double * x, const int incx)
{
	const long double one = 1.0L, zero = 0.0L;

	long double absxi, norm, scale, ssq;
	int ix;

	if (n < 1 || incx < 1) {
		norm = zero;
	} else if (n == 1) {
		norm = mc_fabsl(mc_blas_vector_at(x, 1));
	} else {
		scale = zero;
		ssq   = one;
		for (ix = 1; incx < 0 ? ix >= (n - 1) * incx + 1 : ix <= (n - 1) * incx + 1; ix += incx) {
			absxi = mc_blas_vector_at(x, ix);
			if (absxi != zero) {
				absxi = mc_fabsl(absxi);
				if (scale < absxi) {
					ssq   = one + (ssq * mc_raise2l(scale / absxi));
					scale = absxi;
				} else {
					ssq = ssq + mc_raise2l(absxi / scale);
				}
			}
		}
		norm = scale * mc_sqrtl(ssq);
	}
	return norm;
}

/* \name
 *    ?nrm2 returns the euclidean norm of a vector.
 *
 * \synopsis
 *    real-floating ?nrm2(n, x, incx)
 *    int     incx, n
 *    complex x(*)
 *
 * \purpose
 *    ?nrm2 returns the euclidean norm of a vector: norm2=sqrt(x_*x).
 *
 * \parameters
 *    [in] n     - int. Specifies the number of elements in the input vector `x`.
 *    [in] x     - complex array of size at least (1+(n-1)*abs(incx)).
 *    [in] incx  - int. Specifies the increment for the elements of `x`, incx must not be zero.
 *
 * \examples
 *
 * \level 1 blas routine.
 *     \author Univ. of Tennessee.
 *     \author Univ. of California Berkeley.
 *     \author Univ. of Colorado Denver.
 *     \author NAG Ltd.
 *     \author Sven Hammarling, Nag Ltd.
 */

#pragma mark - mc_blas_scnrm2 -

MC_TARGET_FUNC float mc_blas_scnrm2(const int n, const mc_complex_float_t * x, const int incx)
{
	const float one = 1.0f, zero = 0.0f;

	float cr, ci, norm, scale, ssq, temp;
	int ix;

	if (n < 1 || incx < 1) {
		norm = zero;
	} else {
		scale = zero;
		ssq   = one;
		for (ix = 1; incx < 0 ? ix >= (n - 1) * incx + 1 : ix <= (n - 1) * incx + 1; ix += incx) {
			cr = mc_cmplxrf(mc_blas_vector_at(x, ix));
			ci = mc_cmplxif(mc_blas_vector_at(x, ix));
			if (cr != zero) {
				temp = mc_fabsf(cr);
				if (scale < temp) {
					ssq   = one + (ssq * mc_raise2f(scale / temp));
					scale = temp;
				} else {
					ssq = ssq + mc_raise2f(temp / scale);
				}
			}
			if (ci != zero) {
				temp = mc_fabsf(ci);
				if (scale < temp) {
					ssq   = one + (ssq * mc_raise2f(scale / temp));
					scale = temp;
				} else {
					ssq = ssq + mc_raise2f(temp / scale);
				}
			}
		}
		norm = scale * mc_sqrtf(ssq);
	}
	return norm;
}

#pragma mark - mc_blas_dznrm2 -

MC_TARGET_FUNC double mc_blas_dznrm2(const int n, const mc_complex_double_t * x, const int incx)
{
	const double one = 1.0, zero = 0.0;

	double zr, zi, norm, scale, ssq, temp;
	int ix;

	if (n < 1 || incx < 1) {
		norm = zero;
	} else {
		scale = zero;
		ssq   = one;
		for (ix = 1; incx < 0 ? ix >= (n - 1) * incx + 1 : ix <= (n - 1) * incx + 1; ix += incx) {
			zr = mc_cmplxr(mc_blas_vector_at(x, ix));
			zi = mc_cmplxi(mc_blas_vector_at(x, ix));
			if (zr != zero) {
				temp = mc_fabs(zr);
				if (scale < temp) {
					ssq   = one + (ssq * mc_raise2(scale / temp));
					scale = temp;
				} else {
					ssq = ssq + mc_raise2(temp / scale);
				}
			}
			if (zi != zero) {
				temp = mc_fabs(zi);
				if (scale < temp) {
					ssq   = one + (ssq * mc_raise2(scale / temp));
					scale = temp;
				} else {
					ssq = ssq + mc_raise2(temp / scale);
				}
			}
		}
		norm = scale * mc_sqrt(ssq);
	}
	return norm;
}

#pragma mark - mc_blas_lqnrm2 -

MC_TARGET_FUNC long double mc_blas_lqnrm2(const int n, const mc_complex_long_double_t * x, const int incx)
{
	const long double one = 1.0L, zero = 0.0L;

	long double qr, qi, norm, scale, ssq, temp;
	int ix;

	if (n < 1 || incx < 1) {
		norm = zero;
	} else {
		scale = zero;
		ssq   = one;
		for (ix = 1; incx < 0 ? ix >= (n - 1) * incx + 1 : ix <= (n - 1) * incx + 1; ix += incx) {
			qr = mc_cmplxrl(mc_blas_vector_at(x, ix));
			qi = mc_cmplxil(mc_blas_vector_at(x, ix));
			if (qr != zero) {
				temp = mc_fabsl(qr);
				if (scale < temp) {
					ssq   = one + (ssq * mc_raise2l(scale / temp));
					scale = temp;
				} else {
					ssq = ssq + mc_raise2l(temp / scale);
				}
			}
			if (qi != zero) {
				temp = mc_fabsl(qi);
				if (scale < temp) {
					ssq   = one + (ssq * mc_raise2l(scale / temp));
					scale = temp;
				} else {
					ssq = ssq + mc_raise2l(temp / scale);
				}
			}
		}
		norm = scale * mc_sqrtl(ssq);
	}
	return norm;
}

#endif /* !MC_BLAS_NRM2_H */

/* EOF */
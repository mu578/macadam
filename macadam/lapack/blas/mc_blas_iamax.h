//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_iamax.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

/* \name
 *    i?amax returns the index of the first element having maximum absolute value.
 *
 * \synopsis
 *    int i?amax(n, x, incx)
 *    int           incx, n
 *    real-floating x(*)
 *
 * \purpose
 *    ?nrm2 returns the index of the first element having maximum absolute value.
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

#include <macadam/lapack/blas/mc_blas_abs1.h>
#include <macadam/details/math/mc_fabs.h>

#ifndef MC_BLAS_IAMAX_H
#define MC_BLAS_IAMAX_H

#pragma mark - mc_blas_isamax -

MC_TARGET_FUNC int mc_blas_isamax(const int n, const float * x, const int incx)
{
	float max;
	int i, iamax, ix;

	iamax = 0;
	if (n < 1 || incx <= 0) {
		return iamax - 1;
	}
	iamax = 1;
	if (n == 1) {
		return iamax - 1;
	}
	if (incx == 1) {
		max = mc_fabsf(mc_blas_vector_at(x, 1));
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 2; i <= n; ++i) {
			if (mc_fabsf(mc_blas_vector_at(x, i)) > max) {
				iamax = i;
				max   = mc_fabsf(mc_blas_vector_at(x, i));
			}
		}
	} else {
		ix  = 1;
		max = mc_fabsf(mc_blas_vector_at(x, 1));
		ix  = ix + incx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 2; i <= n; ++i) {
			if (mc_fabsf(mc_blas_vector_at(x, ix)) > max) {
				iamax = i;
				max   = mc_fabsf(mc_blas_vector_at(x, ix));
			}
			ix = ix + incx;
		}
	}
	return iamax - 1;
}

#pragma mark - mc_blas_idamax -

MC_TARGET_FUNC int mc_blas_idamax(const int n, const double * x, const int incx)
{
	double max;
	int i, iamax, ix;

	iamax = 0;
	if (n < 1 || incx <= 0) {
		return iamax - 1;
	}
	iamax = 1;
	if (n == 1) {
		return iamax - 1;
	}
	if (incx == 1) {
		max = mc_fabs(mc_blas_vector_at(x, 1));
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 2; i <= n; ++i) {
			if (mc_fabs(mc_blas_vector_at(x, i)) > max) {
				iamax = i;
				max   = mc_fabs(mc_blas_vector_at(x, i));
			}
		}
	} else {
		ix  = 1;
		max = mc_fabs(mc_blas_vector_at(x, 1));
		ix  = ix + incx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 2; i <= n; ++i) {
			if (mc_fabs(mc_blas_vector_at(x, ix)) > max) {
				iamax = i;
				max   = mc_fabs(mc_blas_vector_at(x, ix));
			}
			ix = ix + incx;
		}
	}
	return iamax - 1;
}

#pragma mark - mc_blas_ilamax -

MC_TARGET_FUNC int mc_blas_ilamax(const int n, const long double * x, const int incx)
{
	long double max;
	int i, iamax, ix;

	iamax = 0;
	if (n < 1 || incx <= 0) {
		return iamax - 1;
	}
	iamax = 1;
	if (n == 1) {
		return iamax - 1;
	}
	if (incx == 1) {
		max = mc_fabsl(mc_blas_vector_at(x, 1));
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 2; i <= n; ++i) {
			if (mc_fabsl(mc_blas_vector_at(x, i)) > max) {
				iamax = i;
				max   = mc_fabsl(mc_blas_vector_at(x, i));
			}
		}
	} else {
		ix  = 1;
		max = mc_fabsl(mc_blas_vector_at(x, 1));
		ix  = ix + incx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 2; i <= n; ++i) {
			if (mc_fabsl(mc_blas_vector_at(x, ix)) > max) {
				iamax = i;
				max   = mc_fabsl(mc_blas_vector_at(x, ix));
			}
			ix = ix + incx;
		}
	}
	return iamax - 1;
}

/* \name
 *    i?amax returns the index of the first element having maximum absolute value.
 *
 * \synopsis
 *    int i?amax(n, x, incx)
 *    int     incx, n
 *    complex x(*)
 *
 * \purpose
 *    ?nrm2 returns the index of the first element having maximum absolute value.
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
 */

#pragma mark - mc_blas_icamax -

MC_TARGET_FUNC int mc_blas_icamax(const int n, const mc_complex_float_t * x, const int incx)
{
	float max;
	int i, iamax, ix;

	iamax = 0;
	if (n < 1 || incx <= 0) {
		return iamax - 1;
	}
	iamax = 1;
	if (n == 1) {
		return iamax - 1;
	}
	if (incx == 1) {
		max = mc_blas_scabs1(mc_blas_vector_at(x, 1));
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 2; i <= n; ++i) {
			if (mc_blas_scabs1(mc_blas_vector_at(x, i)) > max) {
				iamax = i;
				max   = mc_blas_scabs1(mc_blas_vector_at(x, i));
			}
		}
	} else {
		ix  = 1;
		max = mc_blas_scabs1(mc_blas_vector_at(x, 1));
		ix  = ix + incx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 2; i <= n; ++i) {
			if (mc_blas_scabs1(mc_blas_vector_at(x, ix)) > max) {
				iamax = i;
				max   = mc_blas_scabs1(mc_blas_vector_at(x, ix));
			}
			ix = ix + incx;
		}
	}
	return iamax - 1;
}

#pragma mark - mc_blas_izamax -

MC_TARGET_FUNC int mc_blas_izamax(const int n, const mc_complex_double_t * x, const int incx)
{
	double max;
	int i, iamax, ix;

	iamax = 0;
	if (n < 1 || incx <= 0) {
		return iamax - 1;
	}
	iamax = 1;
	if (n == 1) {
		return iamax - 1;
	}
	if (incx == 1) {
		max = mc_blas_dzabs1(mc_blas_vector_at(x, 1));
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 2; i <= n; ++i) {
			if (mc_blas_dzabs1(mc_blas_vector_at(x, i)) > max) {
				iamax = i;
				max   = mc_blas_dzabs1(mc_blas_vector_at(x, i));
			}
		}
	} else {
		ix  = 1;
		max = mc_blas_dzabs1(mc_blas_vector_at(x, 1));
		ix  = ix + incx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 2; i <= n; ++i) {
			if (mc_blas_dzabs1(mc_blas_vector_at(x, ix)) > max) {
				iamax = i;
				max   = mc_blas_dzabs1(mc_blas_vector_at(x, ix));
			}
			ix = ix + incx;
		}
	}
	return iamax - 1;
}

#pragma mark - mc_blas_iqamax -

MC_TARGET_FUNC int mc_blas_iqamax(const int n, const mc_complex_long_double_t * x, const int incx)
{
	long double max;
	int i, iamax, ix;

	iamax = 0;
	if (n < 1 || incx <= 0) {
		return iamax - 1;
	}
	iamax = 1;
	if (n == 1) {
		return iamax - 1;
	}
	if (incx == 1) {
		max = mc_blas_lqabs1(mc_blas_vector_at(x, 1));
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 2; i <= n; ++i) {
			if (mc_blas_lqabs1(mc_blas_vector_at(x, i)) > max) {
				iamax = i;
				max   = mc_blas_lqabs1(mc_blas_vector_at(x, i));
			}
		}
	} else {
		ix  = 1;
		max = mc_blas_lqabs1(mc_blas_vector_at(x, 1));
		ix  = ix + incx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
		for (i = 2; i <= n; ++i) {
			if (mc_blas_lqabs1(mc_blas_vector_at(x, ix)) > max) {
				iamax = i;
				max   = mc_blas_lqabs1(mc_blas_vector_at(x, ix));
			}
			ix = ix + incx;
		}
	}
	return iamax - 1;
}

#endif /* !MC_BLAS_IAMAX_H */

/* EOF */
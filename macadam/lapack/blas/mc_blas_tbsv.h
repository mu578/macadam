//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_tbsv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/lapack/blas/mc_blas_lsame.h>
#include <macadam/lapack/blas/mc_blas_xerbla.h>
#include <macadam/details/math/mc_maxmag.h>
#include <macadam/details/math/mc_minmag.h>

#ifndef MC_BLAS_TBSV_H
#define MC_BLAS_TBSV_H

#pragma mark - mc_blas_stbsv -

MC_TARGET_FUNC void mc_blas_stbsv(const char uplo, const char trans, const char diag, const int n, const int k, const float * a, const int lda, float * x, const int incx)
{
	const float zero = 0.0f;

	float temp;
	int i, info, ix, j, jx, kplus1, kx, l;
	int nounit;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T') && !mc_blas_lsame(trans, 'C')) {
		info = 2;
	} else if (!mc_blas_lsame(diag, 'U') && !mc_blas_lsame(diag, 'N')) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (k < 0) {
		info = 5;
	} else if (lda < k + 1) {
		info = 7;
	} else if (incx == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("STBSV ", info);
		return;
	}

	if (n == 0) {
		return;
	}
	nounit = mc_blas_lsame(diag, 'N');

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					if (mc_blas_vector_at(x, j) != zero) {
						l = kplus1 - j;
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_matrix_at(a, lda, n, kplus1, j);
						}
						temp = mc_blas_vector_at(x, j);
						for (i = (j - 1); i >= mc_maxmag(1, j - k); --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
						}
					}
				}
			} else {
				kx = kx + ((n - 1) * incx);
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					kx = kx - incx;
					if (mc_blas_vector_at(x, jx) != zero) {
						ix = kx;
						l  = kplus1 - j;
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_matrix_at(a, lda, n, kplus1, j);
						}
						temp = mc_blas_vector_at(x, jx);
						for (i = (j - 1); i >= mc_maxmag(1, j - k); --i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
							ix                       = ix - incx;
						}
					}
					jx = jx - incx;
				}
			}
		} else {
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					if (mc_blas_vector_at(x, j) != zero) {
						l = 1 - j;
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_matrix_at(a, lda, n, 1, j);
						}
						temp = mc_blas_vector_at(x, j);
						for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
						}
					}
				}
			} else {
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					kx = kx + incx;
					if (mc_blas_vector_at(x, jx) != zero) {
						ix = kx;
						l  = 1 - j;
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_matrix_at(a, lda, n, 1, j);
						}
						temp = mc_blas_vector_at(x, jx);
						for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
							ix                       = ix + incx;
						}
					}
					jx = jx + incx;
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					temp = mc_blas_vector_at(x, j);
					l    = kplus1 - j;
					for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, kplus1, j);
					}
					mc_blas_vector_at(x, j) = temp;
				}
			} else {
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					temp = mc_blas_vector_at(x, jx);
					ix   = kx;
					l    = kplus1 - j;
					for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
						ix   = ix + incx;
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, kplus1, j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx = jx + incx;
					if (j > k) {
						kx = kx + incx;
					}
				}
			}
		} else {
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, j);
					l    = 1 - j;
					for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, 1, j);
					}
					mc_blas_vector_at(x, j) = temp;
				}
			} else {
				kx = kx + ((n - 1) * incx);
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, jx);
					ix   = kx;
					l    = 1 - j;
					for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
						ix   = ix - incx;
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, 1, j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx = jx - incx;
					if (n - j >= k) {
						kx = kx - incx;
					}
				}
			}
		}
	}
}

#pragma mark - mc_blas_dtbsv -

MC_TARGET_FUNC void mc_blas_dtbsv(const char uplo, const char trans, const char diag, const int n, const int k, const double * a, const int lda, double * x, const int incx)
{
	const double zero = 0.0;

	double temp;
	int i, info, ix, j, jx, kplus1, kx, l;
	int nounit;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T') && !mc_blas_lsame(trans, 'C')) {
		info = 2;
	} else if (!mc_blas_lsame(diag, 'U') && !mc_blas_lsame(diag, 'N')) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (k < 0) {
		info = 5;
	} else if (lda < k + 1) {
		info = 7;
	} else if (incx == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("DTBSV ", info);
		return;
	}

	if (n == 0) {
		return;
	}
	nounit = mc_blas_lsame(diag, 'N');

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					if (mc_blas_vector_at(x, j) != zero) {
						l = kplus1 - j;
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_matrix_at(a, lda, n, kplus1, j);
						}
						temp = mc_blas_vector_at(x, j);
						for (i = (j - 1); i >= mc_maxmag(1, j - k); --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
						}
					}
				}
			} else {
				kx = kx + ((n - 1) * incx);
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					kx = kx - incx;
					if (mc_blas_vector_at(x, jx) != zero) {
						ix = kx;
						l  = kplus1 - j;
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_matrix_at(a, lda, n, kplus1, j);
						}
						temp = mc_blas_vector_at(x, jx);
						for (i = (j - 1); i >= mc_maxmag(1, j - k); --i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
							ix                       = ix - incx;
						}
					}
					jx = jx - incx;
				}
			}
		} else {
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					if (mc_blas_vector_at(x, j) != zero) {
						l = 1 - j;
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_matrix_at(a, lda, n, 1, j);
						}
						temp = mc_blas_vector_at(x, j);
						for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
						}
					}
				}
			} else {
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					kx = kx + incx;
					if (mc_blas_vector_at(x, jx) != zero) {
						ix = kx;
						l  = 1 - j;
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_matrix_at(a, lda, n, 1, j);
						}
						temp = mc_blas_vector_at(x, jx);
						for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
							ix                       = ix + incx;
						}
					}
					jx = jx + incx;
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					temp = mc_blas_vector_at(x, j);
					l    = kplus1 - j;
					for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, kplus1, j);
					}
					mc_blas_vector_at(x, j) = temp;
				}
			} else {
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					temp = mc_blas_vector_at(x, jx);
					ix   = kx;
					l    = kplus1 - j;
					for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
						ix   = ix + incx;
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, kplus1, j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx = jx + incx;
					if (j > k) {
						kx = kx + incx;
					}
				}
			}
		} else {
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, j);
					l    = 1 - j;
					for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, 1, j);
					}
					mc_blas_vector_at(x, j) = temp;
				}
			} else {
				kx = kx + ((n - 1) * incx);
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, jx);
					ix   = kx;
					l    = 1 - j;
					for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
						ix   = ix - incx;
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, 1, j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx = jx - incx;
					if (n - j >= k) {
						kx = kx - incx;
					}
				}
			}
		}
	}
}

#pragma mark - mc_blas_ltbsv -

MC_TARGET_FUNC void mc_blas_ltbsv(const char uplo, const char trans, const char diag, const int n, const int k, const long double * a, const int lda, long double * x, const int incx)
{
	const long double zero = 0.0L;

	long double temp;
	int i, info, ix, j, jx, kplus1, kx, l;
	int nounit;

	info = 0;
	if (!mc_blas_lsame(uplo, 'U') && !mc_blas_lsame(uplo, 'L')) {
		info = 1;
	} else if (!mc_blas_lsame(trans, 'N') && !mc_blas_lsame(trans, 'T') && !mc_blas_lsame(trans, 'C')) {
		info = 2;
	} else if (!mc_blas_lsame(diag, 'U') && !mc_blas_lsame(diag, 'N')) {
		info = 3;
	} else if (n < 0) {
		info = 4;
	} else if (k < 0) {
		info = 5;
	} else if (lda < k + 1) {
		info = 7;
	} else if (incx == 0) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("LTBSV ", info);
		return;
	}

	if (n == 0) {
		return;
	}
	nounit = mc_blas_lsame(diag, 'N');

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	} else if (incx != 1) {
		kx = 1;
	}

	if (mc_blas_lsame(trans, 'N')) {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					if (mc_blas_vector_at(x, j) != zero) {
						l = kplus1 - j;
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_matrix_at(a, lda, n, kplus1, j);
						}
						temp = mc_blas_vector_at(x, j);
						for (i = (j - 1); i >= mc_maxmag(1, j - k); --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
						}
					}
				}
			} else {
				kx = kx + ((n - 1) * incx);
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					kx = kx - incx;
					if (mc_blas_vector_at(x, jx) != zero) {
						ix = kx;
						l  = kplus1 - j;
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_matrix_at(a, lda, n, kplus1, j);
						}
						temp = mc_blas_vector_at(x, jx);
						for (i = (j - 1); i >= mc_maxmag(1, j - k); --i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
							ix                       = ix - incx;
						}
					}
					jx = jx - incx;
				}
			}
		} else {
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					if (mc_blas_vector_at(x, j) != zero) {
						l = 1 - j;
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_matrix_at(a, lda, n, 1, j);
						}
						temp = mc_blas_vector_at(x, j);
						for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
						}
					}
				}
			} else {
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					kx = kx + incx;
					if (mc_blas_vector_at(x, jx) != zero) {
						ix = kx;
						l  = 1 - j;
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_matrix_at(a, lda, n, 1, j);
						}
						temp = mc_blas_vector_at(x, jx);
						for (i = (j + 1); i <= mc_minmag(n, j + k); ++i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_matrix_at(a, lda, n, l + i, j));
							ix                       = ix + incx;
						}
					}
					jx = jx + incx;
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
			kplus1 = k + 1;
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					temp = mc_blas_vector_at(x, j);
					l    = kplus1 - j;
					for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, kplus1, j);
					}
					mc_blas_vector_at(x, j) = temp;
				}
			} else {
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					temp = mc_blas_vector_at(x, jx);
					ix   = kx;
					l    = kplus1 - j;
					for (i = mc_maxmag(1, j - k); i <= (j - 1); ++i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
						ix   = ix + incx;
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, kplus1, j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx = jx + incx;
					if (j > k) {
						kx = kx + incx;
					}
				}
			}
		} else {
			if (incx == 1) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, j);
					l    = 1 - j;
					for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, i));
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, 1, j);
					}
					mc_blas_vector_at(x, j) = temp;
				}
			} else {
				kx = kx + ((n - 1) * incx);
				jx = kx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, jx);
					ix   = kx;
					l    = 1 - j;
					for (i = mc_minmag(n, j + k); i >= (j + 1); --i) {
						temp = temp - (mc_blas_matrix_at(a, lda, n, l + i, j) * mc_blas_vector_at(x, ix));
						ix   = ix - incx;
					}
					if (nounit) {
						temp = temp / mc_blas_matrix_at(a, lda, n, 1, j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx = jx - incx;
					if (n - j >= k) {
						kx = kx - incx;
					}
				}
			}
		}
	}
}

#endif /* !MC_BLAS_TBSV_H */

/* EOF */
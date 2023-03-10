//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_trmv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/lapack/blas/mc_blas_lsame.h>
#include <macadam/lapack/blas/mc_blas_xerbla.h>
#include <macadam/details/math/mc_maxmag.h>

#ifndef MC_BLAS_TRMV_H
#define MC_BLAS_TRMV_H

#pragma mark - mc_blas_strmv -

MC_TARGET_FUNC void mc_blas_strmv(const char uplo, const char trans, const char diag, const int n, const float * a, const int lda, float * x, const int incx)
{
	const float zero = 0.0f;

	float temp;
	int i, info, ix, j, jx, kx;
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
	} else if (lda < mc_maxmag(1, n)) {
		info = 6;
	} else if (incx == 0) {
		info = 8;
	}
	if (info != 0) {
		mc_blas_xerbla("STRMV ", info);
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
						temp = mc_blas_vector_at(x, j);
						for (i = 1; i <= (j - 1); ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) + (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) * mc_blas_matrix_at(a, lda, n, j, j);
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
					if (mc_blas_vector_at(x, jx) != zero) {
						temp = mc_blas_vector_at(x, jx);
						ix   = kx;
						for (i = 1; i <= (j - 1); ++i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) + (temp * mc_blas_matrix_at(a, lda, n, i, j));
							ix                       = ix + incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) * mc_blas_matrix_at(a, lda, n, j, j);
						}
					}
					jx = jx + incx;
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
					if (mc_blas_vector_at(x, j) != zero) {
						temp = mc_blas_vector_at(x, j);
						for (i = n; i >= (j + 1); --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) + (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) * mc_blas_matrix_at(a, lda, n, j, j);
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
					if (mc_blas_vector_at(x, jx) != zero) {
						temp = mc_blas_vector_at(x, jx);
						ix   = kx;
						for (i = n; i >= (j + 1); --i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) + (temp * mc_blas_matrix_at(a, lda, n, i, j));
							ix                       = ix - incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) * mc_blas_matrix_at(a, lda, n, j, j);
						}
					}
					jx = jx + incx;
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
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
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, j, j);
					}
					for (i = (j - 1); i >= 1; --i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, i));
					}
					mc_blas_vector_at(x, j) = temp;
				}
			} else {
				jx = kx + (n - 1) * incx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, jx);
					ix   = jx;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, j, j);
					}
					for (i = (j - 1); i >= 1; --i) {
						ix   = ix + incx;
						temp = temp + (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, ix));
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
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
					temp = mc_blas_vector_at(x, j);
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, j, j);
					}
					for (i = j + 1; i <= n; ++i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, i));
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
					ix   = jx;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, j, j);
					}
					for (i = j + 1; i <= n; ++i) {
						ix   = ix + incx;
						temp = temp + (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, ix));
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
				}
			}
		}
	}
}

#pragma mark - mc_blas_dtrmv -

MC_TARGET_FUNC void mc_blas_dtrmv(const char uplo, const char trans, const char diag, const int n, const double * a, const int lda, double * x, const int incx)
{
	const double zero = 0.0;

	double temp;
	int i, info, ix, j, jx, kx;
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
	} else if (lda < mc_maxmag(1, n)) {
		info = 6;
	} else if (incx == 0) {
		info = 8;
	}
	if (info != 0) {
		mc_blas_xerbla("DTRMV ", info);
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
						temp = mc_blas_vector_at(x, j);
						for (i = 1; i <= (j - 1); ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) + (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) * mc_blas_matrix_at(a, lda, n, j, j);
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
					if (mc_blas_vector_at(x, jx) != zero) {
						temp = mc_blas_vector_at(x, jx);
						ix   = kx;
						for (i = 1; i <= (j - 1); ++i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) + (temp * mc_blas_matrix_at(a, lda, n, i, j));
							ix                       = ix + incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) * mc_blas_matrix_at(a, lda, n, j, j);
						}
					}
					jx = jx + incx;
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
					if (mc_blas_vector_at(x, j) != zero) {
						temp = mc_blas_vector_at(x, j);
						for (i = n; i >= (j + 1); --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) + (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) * mc_blas_matrix_at(a, lda, n, j, j);
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
					if (mc_blas_vector_at(x, jx) != zero) {
						temp = mc_blas_vector_at(x, jx);
						ix   = kx;
						for (i = n; i >= (j + 1); --i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) + (temp * mc_blas_matrix_at(a, lda, n, i, j));
							ix                       = ix - incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) * mc_blas_matrix_at(a, lda, n, j, j);
						}
					}
					jx = jx + incx;
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
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
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, j, j);
					}
					for (i = (j - 1); i >= 1; --i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, i));
					}
					mc_blas_vector_at(x, j) = temp;
				}
			} else {
				jx = kx + (n - 1) * incx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, jx);
					ix   = jx;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, j, j);
					}
					for (i = (j - 1); i >= 1; --i) {
						ix   = ix + incx;
						temp = temp + (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, ix));
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
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
					temp = mc_blas_vector_at(x, j);
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, j, j);
					}
					for (i = j + 1; i <= n; ++i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, i));
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
					ix   = jx;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, j, j);
					}
					for (i = j + 1; i <= n; ++i) {
						ix   = ix + incx;
						temp = temp + (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, ix));
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
				}
			}
		}
	}
}

#pragma mark - mc_blas_ltrmv -

MC_TARGET_FUNC void mc_blas_ltrmv(const char uplo, const char trans, const char diag, const int n, const long double * a, const int lda, long double * x, const int incx)
{
	const long double zero = 0.0L;

	long double temp;
	int i, info, ix, j, jx, kx;
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
	} else if (lda < mc_maxmag(1, n)) {
		info = 6;
	} else if (incx == 0) {
		info = 8;
	}
	if (info != 0) {
		mc_blas_xerbla("LTRMV ", info);
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
						temp = mc_blas_vector_at(x, j);
						for (i = 1; i <= (j - 1); ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) + (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) * mc_blas_matrix_at(a, lda, n, j, j);
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
					if (mc_blas_vector_at(x, jx) != zero) {
						temp = mc_blas_vector_at(x, jx);
						ix   = kx;
						for (i = 1; i <= (j - 1); ++i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) + (temp * mc_blas_matrix_at(a, lda, n, i, j));
							ix                       = ix + incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) * mc_blas_matrix_at(a, lda, n, j, j);
						}
					}
					jx = jx + incx;
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
					if (mc_blas_vector_at(x, j) != zero) {
						temp = mc_blas_vector_at(x, j);
						for (i = n; i >= (j + 1); --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) + (temp * mc_blas_matrix_at(a, lda, n, i, j));
						}
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) * mc_blas_matrix_at(a, lda, n, j, j);
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
					if (mc_blas_vector_at(x, jx) != zero) {
						temp = mc_blas_vector_at(x, jx);
						ix   = kx;
						for (i = n; i >= (j + 1); --i) {
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) + (temp * mc_blas_matrix_at(a, lda, n, i, j));
							ix                       = ix - incx;
						}
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) * mc_blas_matrix_at(a, lda, n, j, j);
						}
					}
					jx = jx + incx;
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
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
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, j, j);
					}
					for (i = (j - 1); i >= 1; --i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, i));
					}
					mc_blas_vector_at(x, j) = temp;
				}
			} else {
				jx = kx + (n - 1) * incx;
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					temp = mc_blas_vector_at(x, jx);
					ix   = jx;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, j, j);
					}
					for (i = (j - 1); i >= 1; --i) {
						ix   = ix + incx;
						temp = temp + (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, ix));
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
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
					temp = mc_blas_vector_at(x, j);
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, j, j);
					}
					for (i = j + 1; i <= n; ++i) {
						temp = temp + (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, i));
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
					ix   = jx;
					if (nounit) {
						temp = temp * mc_blas_matrix_at(a, lda, n, j, j);
					}
					for (i = j + 1; i <= n; ++i) {
						ix   = ix + incx;
						temp = temp + (mc_blas_matrix_at(a, lda, n, i, j) * mc_blas_vector_at(x, ix));
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
				}
			}
		}
	}
}

#endif /* !MC_BLAS_TRMV_H */

/* EOF */
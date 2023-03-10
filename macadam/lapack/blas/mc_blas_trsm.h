//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_trsm.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/lapack/blas/mc_blas_lsame.h>
#include <macadam/lapack/blas/mc_blas_xerbla.h>
#include <macadam/details/math/mc_maxmag.h>

#ifndef MC_BLAS_TRSM_H
#define MC_BLAS_TRSM_H

#pragma mark - mc_blas_strsm -

MC_TARGET_FUNC void mc_blas_strsm(const char side, const char uplo, const char transa, const char diag, const int m, const int n, const float alpha, const float * a, const int lda, float * b, const int ldb)
{
	const float one = 1.0f, zero = 0.0f;

	float temp;
	int i, info, j, k, nrowa, ka;
	int lside, nounit, upper;

	lside = mc_blas_lsame(side, 'L');
	if (lside) {
		ka    = n;
		nrowa = m;
		mc_unused(ka);
	} else {
		ka    = m;
		nrowa = n;
		mc_unused(ka);
	}
	nounit = mc_blas_lsame(diag, 'N');
	upper  = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (! lside && !mc_blas_lsame(side, 'R')) {
		info = 1;
	} else if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 2;
	} else if (!mc_blas_lsame(transa, 'N') && !mc_blas_lsame(transa, 'T') && !mc_blas_lsame(transa, 'C')) {
		info = 3;
	} else if (!mc_blas_lsame(diag, 'U') && !mc_blas_lsame(diag, 'N')) {
		info = 4;
	} else if (m < 0) {
		info = 5;
	} else if (n < 0) {
		info = 6;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 9;
	} else if (ldb < mc_maxmag(1, m)) {
		info = 11;
	}
	if (info != 0) {
		mc_blas_xerbla("STRSM ", info);
		return;
	}

	if (m == 0 || n == 0) {
		return;
	}

	if (alpha == zero) {
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= m; ++i) {
				mc_blas_matrix_at(b, ldb, n, i, j) = zero;
			}
		}
		return;
	}

	if (lside) {
		if (mc_blas_lsame(transa, 'N')) {
			if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
					for (k = m; k >= 1; --k) {
						if (mc_blas_matrix_at(b, ldb, n, k, j) != zero) {
							if (nounit) {
								mc_blas_matrix_at(b, ldb, n, k, j) = mc_blas_matrix_at(b, ldb, n, k, j) / mc_blas_matrix_at(a, lda, ka, k, k);
							}
							for (i = 1; i <= (k - 1); ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (mc_blas_matrix_at(b, ldb, n, k, j) * mc_blas_matrix_at(a, lda, ka, i, k));
							}
						}
					}
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
					for (k = 1; k <= m; ++k) {
						if (mc_blas_matrix_at(b, ldb, n, k, j) != zero) {
							if (nounit) {
								mc_blas_matrix_at(b, ldb, n, k, j) = mc_blas_matrix_at(b, ldb, n, k, j) / mc_blas_matrix_at(a, lda, ka, k, k);
							}
							for (i = (k + 1); i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (mc_blas_matrix_at(b, ldb, n, k, j) * mc_blas_matrix_at(a, lda, ka, i, k));
							}
						}
					}
				}
			}
		} else {
			if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = 1; i <= m; ++i) {
						temp = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						for (k = 1; k <= (i - 1); ++k) {
							temp = temp - (mc_blas_matrix_at(a, lda, ka, k, i) * mc_blas_matrix_at(b, ldb, n, k, j));
						}
						if (nounit) {
							temp = temp / mc_blas_matrix_at(a, lda, ka, i, i);
						}
						mc_blas_matrix_at(b, ldb, n, i, j) = temp;
					}
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = m; i >= 1; --i) {
						temp = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						for (k = i + 1; k <= m; ++k) {
							temp = temp - (mc_blas_matrix_at(a, lda, ka, k, i) * mc_blas_matrix_at(b, ldb, n, k, j));
						}
						if (nounit) {
							temp = temp / mc_blas_matrix_at(a, lda, ka, i, i);
						}
						mc_blas_matrix_at(b, ldb, n, i, j) = temp;
					}
				}
			}
		}
	} else {
		if (mc_blas_lsame(transa, 'N')) {
			if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
					for (k = 1; k <= (j - 1); ++k) {
						if (mc_blas_matrix_at(a, lda, ka, k, j) != zero) {
							for (i = 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (mc_blas_matrix_at(a, lda, ka, k, j) * mc_blas_matrix_at(b, ldb, n, i, k));
							}
						}
					}
					if (nounit) {
						temp = one / mc_blas_matrix_at(a, lda, ka, j, j);
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = temp * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
					for (k = j + 1; k <= n; ++k) {
						if (mc_blas_matrix_at(a, lda, ka, k, j) != zero) {
							for (i = 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (mc_blas_matrix_at(a, lda, ka, k, j) * mc_blas_matrix_at(b, ldb, n, i, k));
							}
						}
					}
					if (nounit) {
						temp = one / mc_blas_matrix_at(a, lda, ka, j, j);
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = temp * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
				}
			}
		} else {
			if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (k = n; k >= 1; --k) {
					if (nounit) {
						temp = one / mc_blas_matrix_at(a, lda, ka, k, k);
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, k) = temp * mc_blas_matrix_at(b, ldb, n, i, k);
						}
					}
					for (j = 1; j <= (k - 1); ++j) {
						if (mc_blas_matrix_at(a, lda, ka, j, k) != zero) {
							temp = mc_blas_matrix_at(a, lda, ka, j, k);
							for (i = 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (temp * mc_blas_matrix_at(b, ldb, n, i, k));
							}
						}
					}
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, k) = alpha * mc_blas_matrix_at(b, ldb, n, i, k);
						}
					}
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (k = 1; k <= n; ++k) {
					if (nounit) {
						temp = one / mc_blas_matrix_at(a, lda, ka, k, k);
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, k) = temp * mc_blas_matrix_at(b, ldb, n, i, k);
						}
					}
					for (j = (k + 1); j <= n; ++j) {
						if (mc_blas_matrix_at(a, lda, ka, j, k) != zero) {
							temp = mc_blas_matrix_at(a, lda, ka, j, k);
							for (i = 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (temp * mc_blas_matrix_at(b, ldb, n, i, k));
							}
						}
					}
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, k) = alpha * mc_blas_matrix_at(b, ldb, n, i, k);
						}
					}
				}
			}
		}
	}
}

#pragma mark - mc_blas_dtrsm -

MC_TARGET_FUNC void mc_blas_dtrsm(const char side, const char uplo, const char transa, const char diag, const int m, const int n, const double alpha, const double * a, const int lda, double * b, const int ldb)
{
	const double one = 1.0, zero = 0.0;

	double temp;
	int i, info, j, k, nrowa, ka;
	int lside, nounit, upper;

	lside = mc_blas_lsame(side, 'L');
	if (lside) {
		ka    = n;
		nrowa = m;
		mc_unused(ka);
	} else {
		ka    = m;
		nrowa = n;
		mc_unused(ka);
	}
	nounit = mc_blas_lsame(diag, 'N');
	upper  = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (! lside && !mc_blas_lsame(side, 'R')) {
		info = 1;
	} else if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 2;
	} else if (!mc_blas_lsame(transa, 'N') && !mc_blas_lsame(transa, 'T') && !mc_blas_lsame(transa, 'C')) {
		info = 3;
	} else if (!mc_blas_lsame(diag, 'U') && !mc_blas_lsame(diag, 'N')) {
		info = 4;
	} else if (m < 0) {
		info = 5;
	} else if (n < 0) {
		info = 6;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 9;
	} else if (ldb < mc_maxmag(1, m)) {
		info = 11;
	}
	if (info != 0) {
		mc_blas_xerbla("DTRSM ", info);
		return;
	}

	if (m == 0 || n == 0) {
		return;
	}

	if (alpha == zero) {
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= m; ++i) {
				mc_blas_matrix_at(b, ldb, n, i, j) = zero;
			}
		}
		return;
	}

	if (lside) {
		if (mc_blas_lsame(transa, 'N')) {
			if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
					for (k = m; k >= 1; --k) {
						if (mc_blas_matrix_at(b, ldb, n, k, j) != zero) {
							if (nounit) {
								mc_blas_matrix_at(b, ldb, n, k, j) = mc_blas_matrix_at(b, ldb, n, k, j) / mc_blas_matrix_at(a, lda, ka, k, k);
							}
							for (i = 1; i <= (k - 1); ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (mc_blas_matrix_at(b, ldb, n, k, j) * mc_blas_matrix_at(a, lda, ka, i, k));
							}
						}
					}
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
					for (k = 1; k <= m; ++k) {
						if (mc_blas_matrix_at(b, ldb, n, k, j) != zero) {
							if (nounit) {
								mc_blas_matrix_at(b, ldb, n, k, j) = mc_blas_matrix_at(b, ldb, n, k, j) / mc_blas_matrix_at(a, lda, ka, k, k);
							}
							for (i = k + 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (mc_blas_matrix_at(b, ldb, n, k, j) * mc_blas_matrix_at(a, lda, ka, i, k));
							}
						}
					}
				}
			}
		} else {
			if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = 1; i <= m; ++i) {
						temp = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						for (k = 1; k <= (i - 1); ++k) {
							temp = temp - (mc_blas_matrix_at(a, lda, ka, k, i) * mc_blas_matrix_at(b, ldb, n, k, j));
						}
						if (nounit) {
							temp = temp / mc_blas_matrix_at(a, lda, ka, i, i);
						}
						mc_blas_matrix_at(b, ldb, n, i, j) = temp;
					}
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = m; i >= 1; --i) {
						temp = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						for (k = i + 1; k <= m; ++k) {
							temp = temp - (mc_blas_matrix_at(a, lda, ka, k, i) * mc_blas_matrix_at(b, ldb, n, k, j));
						}
						if (nounit) {
							temp = temp / mc_blas_matrix_at(a, lda, ka, i, i);
						}
						mc_blas_matrix_at(b, ldb, n, i, j) = temp;
					}
				}
			}
		}
	} else {
		if (mc_blas_lsame(transa, 'N')) {
			if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
					for (k = 1; k <= (j - 1); ++k) {
						if (mc_blas_matrix_at(a, lda, ka, k, j) != zero) {
							for (i = 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (mc_blas_matrix_at(a, lda, ka, k, j) * mc_blas_matrix_at(b, ldb, n, i, k));
							}
						}
					}
					if (nounit) {
						temp = one / mc_blas_matrix_at(a, lda, ka, j, j);
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = temp * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
					for (k = j + 1; k <= n; ++k) {
						if (mc_blas_matrix_at(a, lda, ka, k, j) != zero) {
							for (i = 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (mc_blas_matrix_at(a, lda, ka, k, j) * mc_blas_matrix_at(b, ldb, n, i, k));
							}
						}
					}
					if (nounit) {
						temp = one / mc_blas_matrix_at(a, lda, ka, j, j);
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = temp * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
				}
			}
		} else {
			if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (k = n; k >= 1; --k) {
					if (nounit) {
						temp = one / mc_blas_matrix_at(a, lda, ka, k, k);
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, k) = temp * mc_blas_matrix_at(b, ldb, n, i, k);
						}
					}
					for (j = 1; j <= (k - 1); ++j) {
						if (mc_blas_matrix_at(a, lda, ka, j, k) != zero) {
							temp = mc_blas_matrix_at(a, lda, ka, j, k);
							for (i = 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (temp * mc_blas_matrix_at(b, ldb, n, i, k));
							}
						}
					}
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, k) = alpha * mc_blas_matrix_at(b, ldb, n, i, k);
						}
					}
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (k = 1; k <= n; ++k) {
					if (nounit) {
						temp = one / mc_blas_matrix_at(a, lda, ka, k, k);
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, k) = temp * mc_blas_matrix_at(b, ldb, n, i, k);
						}
					}
					for (j = (k + 1); j <= n; ++j) {
						if (mc_blas_matrix_at(a, lda, ka, j, k) != zero) {
							temp = mc_blas_matrix_at(a, lda, ka, j, k);
							for (i = 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (temp * mc_blas_matrix_at(b, ldb, n, i, k));
							}
						}
					}
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, k) = alpha * mc_blas_matrix_at(b, ldb, n, i, k);
						}
					}
				}
			}
		}
	}
}

#pragma mark - mc_blas_ltrsm -

MC_TARGET_FUNC void mc_blas_ltrsm(const char side, const char uplo, const char transa, const char diag, const int m, const int n, const long double alpha, const long double * a, const int lda, long double * b, const int ldb)
{
	const long double one = 1.0L, zero = 0.0L;

	long double temp;
	int i, info, j, k, nrowa, ka;
	int lside, nounit, upper;

	lside = mc_blas_lsame(side, 'L');
	if (lside) {
		ka    = n;
		nrowa = m;
		mc_unused(ka);
	} else {
		ka    = m;
		nrowa = n;
		mc_unused(ka);
	}
	nounit = mc_blas_lsame(diag, 'N');
	upper  = mc_blas_lsame(uplo, 'U');

	info = 0;
	if (! lside && !mc_blas_lsame(side, 'R')) {
		info = 1;
	} else if (!upper && !mc_blas_lsame(uplo, 'L')) {
		info = 2;
	} else if (!mc_blas_lsame(transa, 'N') && !mc_blas_lsame(transa, 'T') && !mc_blas_lsame(transa, 'C')) {
		info = 3;
	} else if (!mc_blas_lsame(diag, 'U') && !mc_blas_lsame(diag, 'N')) {
		info = 4;
	} else if (m < 0) {
		info = 5;
	} else if (n < 0) {
		info = 6;
	} else if (lda < mc_maxmag(1, nrowa)) {
		info = 9;
	} else if (ldb < mc_maxmag(1, m)) {
		info = 11;
	}
	if (info != 0) {
		mc_blas_xerbla("LTRSM ", info);
		return;
	}

	if (m == 0 || n == 0) {
		return;
	}

	if (alpha == zero) {
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= m; ++i) {
				mc_blas_matrix_at(b, ldb, n, i, j) = zero;
			}
		}
		return;
	}

	if (lside) {
		if (mc_blas_lsame(transa, 'N')) {
			if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
					for (k = m; k >= 1; --k) {
						if (mc_blas_matrix_at(b, ldb, n, k, j) != zero) {
							if (nounit) {
								mc_blas_matrix_at(b, ldb, n, k, j) = mc_blas_matrix_at(b, ldb, n, k, j) / mc_blas_matrix_at(a, lda, ka, k, k);
							}
							for (i = 1; i <= (k - 1); ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (mc_blas_matrix_at(b, ldb, n, k, j) * mc_blas_matrix_at(a, lda, ka, i, k));
							}
						}
					}
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
					for (k = 1; k <= m; ++k) {
						if (mc_blas_matrix_at(b, ldb, n, k, j) != zero) {
							if (nounit) {
								mc_blas_matrix_at(b, ldb, n, k, j) = mc_blas_matrix_at(b, ldb, n, k, j) / mc_blas_matrix_at(a, lda, ka, k, k);
							}
							for (i = k + 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (mc_blas_matrix_at(b, ldb, n, k, j) * mc_blas_matrix_at(a, lda, ka, i, k));
							}
						}
					}
				}
			}
		} else {
			if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = 1; i <= m; ++i) {
						temp = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						for (k = 1; k <= (i - 1); ++k) {
							temp = temp - (mc_blas_matrix_at(a, lda, ka, k, i) * mc_blas_matrix_at(b, ldb, n, k, j));
						}
						if (nounit) {
							temp = temp / mc_blas_matrix_at(a, lda, ka, i, i);
						}
						mc_blas_matrix_at(b, ldb, n, i, j) = temp;
					}
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					for (i = m; i >= 1; --i) {
						temp = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						for (k = i + 1; k <= m; ++k) {
							temp = temp - (mc_blas_matrix_at(a, lda, ka, k, i) * mc_blas_matrix_at(b, ldb, n, k, j));
						}
						if (nounit) {
							temp = temp / mc_blas_matrix_at(a, lda, ka, i, i);
						}
						mc_blas_matrix_at(b, ldb, n, i, j) = temp;
					}
				}
			}
		}
	} else {
		if (mc_blas_lsame(transa, 'N')) {
			if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = 1; j <= n; ++j) {
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
					for (k = 1; k <= (j - 1); ++k) {
						if (mc_blas_matrix_at(a, lda, ka, k, j) != zero) {
							for (i = 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (mc_blas_matrix_at(a, lda, ka, k, j) * mc_blas_matrix_at(b, ldb, n, i, k));
							}
						}
					}
					if (nounit) {
						temp = one / mc_blas_matrix_at(a, lda, ka, j, j);
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = temp * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (j = n; j >= 1; --j) {
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = alpha * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
					for (k = j + 1; k <= n; ++k) {
						if (mc_blas_matrix_at(a, lda, ka, k, j) != zero) {
							for (i = 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (mc_blas_matrix_at(a, lda, ka, k, j) * mc_blas_matrix_at(b, ldb, n, i, k));
							}
						}
					}
					if (nounit) {
						temp = one / mc_blas_matrix_at(a, lda, ka, j, j);
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, j) = temp * mc_blas_matrix_at(b, ldb, n, i, j);
						}
					}
				}
			}
		} else {
			if (upper) {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (k = n; k >= 1; --k) {
					if (nounit) {
						temp = one / mc_blas_matrix_at(a, lda, ka, k, k);
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, k) = temp * mc_blas_matrix_at(b, ldb, n, i, k);
						}
					}
					for (j = 1; j <= (k - 1); ++j) {
						if (mc_blas_matrix_at(a, lda, ka, j, k) != zero) {
							temp = mc_blas_matrix_at(a, lda, ka, j, k);
							for (i = 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (temp * mc_blas_matrix_at(b, ldb, n, i, k));
							}
						}
					}
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, k) = alpha * mc_blas_matrix_at(b, ldb, n, i, k);
						}
					}
				}
			} else {
#	if MC_TARGET_OPENMP
#		if MC_TARGET_OPENMP_LOOP_USE_PARALLEL_FOR
#			pragma omp parallel for
#		elif MC_TARGET_OPENMP_LOOP_USE_FOR_SIMD
#			pragma omp for simd
#		endif
#	endif
				for (k = 1; k <= n; ++k) {
					if (nounit) {
						temp = one / mc_blas_matrix_at(a, lda, ka, k, k);
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, k) = temp * mc_blas_matrix_at(b, ldb, n, i, k);
						}
					}
					for (j = (k + 1); j <= n; ++j) {
						if (mc_blas_matrix_at(a, lda, ka, j, k) != zero) {
							temp = mc_blas_matrix_at(a, lda, ka, j, k);
							for (i = 1; i <= m; ++i) {
								mc_blas_matrix_at(b, ldb, n, i, j) = mc_blas_matrix_at(b, ldb, n, i, j) - (temp * mc_blas_matrix_at(b, ldb, n, i, k));
							}
						}
					}
					if (alpha != one) {
						for (i = 1; i <= m; ++i) {
							mc_blas_matrix_at(b, ldb, n, i, k) = alpha * mc_blas_matrix_at(b, ldb, n, i, k);
						}
					}
				}
			}
		}
	}
}

#endif /* !MC_BLAS_TRSM_H */

/* EOF */
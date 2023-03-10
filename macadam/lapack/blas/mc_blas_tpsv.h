//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_tpsv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/blas/mc_blas_access.h>
#include <macadam/lapack/blas/mc_blas_lsame.h>
#include <macadam/lapack/blas/mc_blas_xerbla.h>

#ifndef MC_BLAS_TPSV_H
#define MC_BLAS_TPSV_H

#pragma mark - mc_blas_stpsv -

MC_TARGET_FUNC void mc_blas_stpsv(const char uplo, const char trans, const char diag, const int n, const float * ap, float * x, const int incx)
{
	const float zero = 0.0f;

	float temp;
	int i, info, ix, j, jx, k, kk, kx;
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
	} else if (incx == 0) {
		info = 7;
	}
	if (info != 0) {
		mc_blas_xerbla("STPSV ", info);
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
			kk = n * (n + 1) / 2;
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
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_vector_at(ap, kk);
						}
						temp = mc_blas_vector_at(x, j);
						k    = kk - 1;
						for (i = (j - 1); i >= 1; --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_vector_at(ap, k));
							k                       = k - 1;
						}
					}
					kk = kk - j;
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
					if (mc_blas_vector_at(x, jx) != zero) {
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_vector_at(ap, kk);
						}
						temp = mc_blas_vector_at(x, jx);
						ix   = jx;
						for (k = (kk - 1); k >= (kk - j + 1); --k) {
							ix                       = ix - incx;
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_vector_at(ap, k));
						}
					}
					jx = jx - incx;
					kk = kk - j;
				}
			}
		} else {
			kk = 1;
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
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_vector_at(ap, kk);
						}
						temp = mc_blas_vector_at(x, j);
						k    = kk + 1;
						for (i = (j + 1); i <= n; ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_vector_at(ap, k));
							k                       = k + 1;
						}
					}
					kk = kk + (n - j + 1);
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
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_vector_at(ap, kk);
						}
						temp = mc_blas_vector_at(x, jx);
						ix = jx;
						for (k = (kk + 1); k <= (kk + n - j); ++k) {
							ix = ix + incx;
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_vector_at(ap, k));
						}
					}
					jx = jx + incx;
					kk = kk + (n - j + 1);
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
			kk = 1;
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
					k    = kk;
					for (i = 1; i <= (j - 1); ++i) {
						temp = temp - (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, i));
						k    = k + 1;
					}
					if (nounit) {
						temp = temp / mc_blas_vector_at(ap, kk + j - 1);
					}
					mc_blas_vector_at(x, j) = temp;
					kk                      = kk + j;
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
					for (k = kk; k <= (kk + j - 2); ++k) {
						temp = temp - (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, ix));
						ix   = ix + incx;
					}
					if (nounit) {
						temp = temp / mc_blas_vector_at(ap, kk + j - 1);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
					kk                       = kk + j;
				}
			}
		} else {
			kk = n * (n + 1) / 2;
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
					k    = kk;
					for (i = n; i >= (j + 1); --i) {
						temp = temp - (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, i));
						k    = k - 1;
					}
					if (nounit) {
						temp = temp / mc_blas_vector_at(ap, kk - n + j);
					}
					mc_blas_vector_at(x, j) = temp;
					kk                      = kk - (n - j + 1);
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
					for (k = kk; k >= (kk - (n - (j + 1))); --k) {
						temp = temp - (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, ix));
						ix   = ix - incx;
					}
					if (nounit) {
						temp = temp / mc_blas_vector_at(ap, kk - n + j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
					kk                       = kk - (n - j + 1);
				}
			}
		}
	}
}

#pragma mark - mc_blas_dtpsv -

MC_TARGET_FUNC void mc_blas_dtpsv(const char uplo, const char trans, const char diag, const int n, const double * ap, double * x, const int incx)
{
	const double zero = 0.0;

	double temp;
	int i, info, ix, j, jx, k, kk, kx;
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
	} else if (incx == 0) {
		info = 7;
	}
	if (info != 0) {
		mc_blas_xerbla("DTPSV ", info);
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
			kk = n * (n + 1) / 2;
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
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_vector_at(ap, kk);
						}
						temp = mc_blas_vector_at(x, j);
						k    = kk - 1;
						for (i = (j - 1); i >= 1; --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_vector_at(ap, k));
							k                       = k - 1;
						}
					}
					kk = kk - j;
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
					if (mc_blas_vector_at(x, jx) != zero) {
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_vector_at(ap, kk);
						}
						temp = mc_blas_vector_at(x, jx);
						ix   = jx;
						for (k = (kk - 1); k >= (kk - j + 1); --k) {
							ix                       = ix - incx;
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_vector_at(ap, k));
						}
					}
					jx = jx - incx;
					kk = kk - j;
				}
			}
		} else {
			kk = 1;
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
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_vector_at(ap, kk);
						}
						temp = mc_blas_vector_at(x, j);
						k    = kk + 1;
						for (i = (j + 1); i <= n; ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_vector_at(ap, k));
							k                       = k + 1;
						}
					}
					kk = kk + (n - j + 1);
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
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_vector_at(ap, kk);
						}
						temp = mc_blas_vector_at(x, jx);
						ix = jx;
						for (k = (kk + 1); k <= (kk + n - j); ++k) {
							ix = ix + incx;
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_vector_at(ap, k));
						}
					}
					jx = jx + incx;
					kk = kk + (n - j + 1);
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
			kk = 1;
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
					k    = kk;
					for (i = 1; i <= (j - 1); ++i) {
						temp = temp - (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, i));
						k    = k + 1;
					}
					if (nounit) {
						temp = temp / mc_blas_vector_at(ap, kk + j - 1);
					}
					mc_blas_vector_at(x, j) = temp;
					kk                      = kk + j;
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
					for (k = kk; k <= (kk + j - 2); ++k) {
						temp = temp - (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, ix));
						ix   = ix + incx;
					}
					if (nounit) {
						temp = temp / mc_blas_vector_at(ap, kk + j - 1);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
					kk                       = kk + j;
				}
			}
		} else {
			kk = n * (n + 1) / 2;
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
					k    = kk;
					for (i = n; i >= (j + 1); --i) {
						temp = temp - (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, i));
						k    = k - 1;
					}
					if (nounit) {
						temp = temp / mc_blas_vector_at(ap, kk - n + j);
					}
					mc_blas_vector_at(x, j) = temp;
					kk                      = kk - (n - j + 1);
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
					for (k = kk; k >= (kk - (n - (j + 1))); --k) {
						temp = temp - (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, ix));
						ix   = ix - incx;
					}
					if (nounit) {
						temp = temp / mc_blas_vector_at(ap, kk - n + j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
					kk                       = kk - (n - j + 1);
				}
			}
		}
	}
}

#pragma mark - mc_blas_ltpsv -

MC_TARGET_FUNC void mc_blas_ltpsv(const char uplo, const char trans, const char diag, const int n, const long double * ap, long double * x, const int incx)
{
	const long double zero = 0.0L;

	long double temp;
	int i, info, ix, j, jx, k, kk, kx;
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
	} else if (incx == 0) {
		info = 7;
	}
	if (info != 0) {
		mc_blas_xerbla("LTPSV ", info);
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
			kk = n * (n + 1) / 2;
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
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_vector_at(ap, kk);
						}
						temp = mc_blas_vector_at(x, j);
						k    = kk - 1;
						for (i = (j - 1); i >= 1; --i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_vector_at(ap, k));
							k                       = k - 1;
						}
					}
					kk = kk - j;
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
					if (mc_blas_vector_at(x, jx) != zero) {
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_vector_at(ap, kk);
						}
						temp = mc_blas_vector_at(x, jx);
						ix   = jx;
						for (k = (kk - 1); k >= (kk - j + 1); --k) {
							ix                       = ix - incx;
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_vector_at(ap, k));
						}
					}
					jx = jx - incx;
					kk = kk - j;
				}
			}
		} else {
			kk = 1;
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
						if (nounit) {
							mc_blas_vector_at(x, j) = mc_blas_vector_at(x, j) / mc_blas_vector_at(ap, kk);
						}
						temp = mc_blas_vector_at(x, j);
						k    = kk + 1;
						for (i = (j + 1); i <= n; ++i) {
							mc_blas_vector_at(x, i) = mc_blas_vector_at(x, i) - (temp * mc_blas_vector_at(ap, k));
							k                       = k + 1;
						}
					}
					kk = kk + (n - j + 1);
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
						if (nounit) {
							mc_blas_vector_at(x, jx) = mc_blas_vector_at(x, jx) / mc_blas_vector_at(ap, kk);
						}
						temp = mc_blas_vector_at(x, jx);
						ix = jx;
						for (k = (kk + 1); k <= (kk + n - j); ++k) {
							ix = ix + incx;
							mc_blas_vector_at(x, ix) = mc_blas_vector_at(x, ix) - (temp * mc_blas_vector_at(ap, k));
						}
					}
					jx = jx + incx;
					kk = kk + (n - j + 1);
				}
			}
		}
	} else {
		if (mc_blas_lsame(uplo, 'U')) {
			kk = 1;
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
					k    = kk;
					for (i = 1; i <= (j - 1); ++i) {
						temp = temp - (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, i));
						k    = k + 1;
					}
					if (nounit) {
						temp = temp / mc_blas_vector_at(ap, kk + j - 1);
					}
					mc_blas_vector_at(x, j) = temp;
					kk                      = kk + j;
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
					for (k = kk; k <= (kk + j - 2); ++k) {
						temp = temp - (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, ix));
						ix   = ix + incx;
					}
					if (nounit) {
						temp = temp / mc_blas_vector_at(ap, kk + j - 1);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx + incx;
					kk                       = kk + j;
				}
			}
		} else {
			kk = n * (n + 1) / 2;
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
					k    = kk;
					for (i = n; i >= (j + 1); --i) {
						temp = temp - (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, i));
						k    = k - 1;
					}
					if (nounit) {
						temp = temp / mc_blas_vector_at(ap, kk - n + j);
					}
					mc_blas_vector_at(x, j) = temp;
					kk                      = kk - (n - j + 1);
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
					for (k = kk; k >= (kk - (n - (j + 1))); --k) {
						temp = temp - (mc_blas_vector_at(ap, k) * mc_blas_vector_at(x, ix));
						ix   = ix - incx;
					}
					if (nounit) {
						temp = temp / mc_blas_vector_at(ap, kk - n + j);
					}
					mc_blas_vector_at(x, jx) = temp;
					jx                       = jx - incx;
					kk                       = kk - (n - j + 1);
				}
			}
		}
	}
}

#endif /* !MC_BLAS_TPSV_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lascl.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_lapack_lamch.h>
#include <macadam/lapack/mc_lapack_lartgp.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_isnan.h>

#ifndef MC_LAPACKE_LASCL_H
#define MC_LAPACKE_LASCL_H

#pragma mark - mc_lapack_slascl -

MC_TARGET_FUNC void mc_lapack_slascl(const char type, const int kl, const int ku, float cfrom, float cto, const int m, const int n, float * a, const int lda, int * info)
{
	const float one = 1.0f, zero = 0.0f;

	int done;
	int i, itype, j, k1, k2, k3, k4;
	float bignum, cfrom1, cfromc, cto1, ctoc, mul, smlnum;

	*info = 0;

	if (mc_blas_lsame(type, 'G')) {
		itype = 0;
	} else if (mc_blas_lsame(type, 'L')) {
		itype = 1;
	} else if (mc_blas_lsame(type, 'U')) {
		itype = 2;
	} else if (mc_blas_lsame(type, 'H')) {
		itype = 3;
	} else if (mc_blas_lsame(type, 'B')) {
		itype = 4;
	} else if (mc_blas_lsame(type, 'Q')) {
		itype = 5;
	} else if (mc_blas_lsame(type, 'Z')) {
		itype = 6;
	} else {
		itype = -1;
	}

	if (itype == -1) {
		*info = -1;
	} else if (cfrom == zero || mc_isnan(cfrom)) {
		*info = -4;
	} else if (mc_isnan(cto)) {
		*info = -5;
	} else if (m < 0) {
		*info = -6;
	} else if (n < 0 || (itype == 4 && n != m) || (itype == 5 && n != m)) {
		*info = -7;
	} else if (itype <= 3 && lda < mc_maxmag(1, m)) {
		*info = -9;
	} else if (itype >= 4) {
		if (kl < 0 || kl > mc_maxmag(m - 1, 0)) {
			*info = -2;
		} else {
			if (
				   (ku < 0 || ku > mc_maxmag(n - 1, 0))
				|| ((itype == 4 || itype == 5) && kl != ku)
			) {
				*info = -3;
			} else if (
				   (itype == 4 && lda < (kl + 1))
				|| (itype == 5 && lda < (ku + 1))
				|| (itype == 6 && lda < (2 * kl + ku + 1))
			) {
				*info = -9;
			}
		}
	}

	if (*info != 0) {
		mc_blas_xerbla("SLASCL", -(*info));
		return;
	}

	if (n == 0 || m == 0) {
		return;
	}

	smlnum = mc_lapack_slamch('S');
	bignum = one / smlnum;

	cfromc = cfrom;
	ctoc   = cto;

F10:
	cfrom1 = cfromc * smlnum;
	if (cfrom1 == cfromc) {
		mul  = ctoc / cfromc;
		done = 1;
		cto1 = ctoc;
		mc_unused(cto1);
	} else {
		cto1 = ctoc / bignum;
		if (cto1 == ctoc) {
			mul    = ctoc;
			done   = 1;
			cfromc = one;
		} else if (mc_fabsf(cfrom1) > mc_fabsf(ctoc) && ctoc != zero) {
			mul    = smlnum;
			done   = 0;
			cfromc = cfrom1;
		} else if (mc_fabsf(cto1) > mc_fabsf(cfromc)) {
			mul  = bignum;
			done = 0;
			ctoc = cto1;
		} else {
			mul  = ctoc / cfromc;
			done = 1;
		}
	}

	if (itype == 0) {
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= m; ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 1) {
		for (j = 1; j <= n; ++j) {
			for (i = j; i <= m; ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 2) {
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= mc_minmag(j, m); ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 3) {
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= mc_minmag(j + 1, m); ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 4) {
		k3 = kl + 1;
		k4 = n + 1;
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= mc_minmag(k3, k4 - j); ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 5) {
		k1 = ku + 2;
		k3 = ku + 1;
		for (j = 1; j <= n; ++j) {
			for (i = mc_maxmag(k1 - j, 1); i <= k3; ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 6) {
		k1 = kl + ku + 2;
		k2 = kl + 1;
		k3 = (2 * kl) + ku + 1;
		k4 = kl + ku + 1 + m;
		for (j = 1; j <= n; ++j) {
			for (i = mc_maxmag(k1 - j, k2); i <=  mc_minmag(k3, k4 - j); ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	}
	if (!done) {
		goto F10;
	}
}

#pragma mark - mc_lapack_dlascl -

MC_TARGET_FUNC void mc_lapack_dlascl(const char type, const int kl, const int ku, double cfrom, double cto, const int m, const int n, double * a, const int lda, int * info)
{
	const double one = 1.0, zero = 0.0;

	int done;
	int i, itype, j, k1, k2, k3, k4;
	double bignum, cfrom1, cfromc, cto1, ctoc, mul, smlnum;

	*info = 0;

	if (mc_blas_lsame(type, 'G')) {
		itype = 0;
	} else if (mc_blas_lsame(type, 'L')) {
		itype = 1;
	} else if (mc_blas_lsame(type, 'U')) {
		itype = 2;
	} else if (mc_blas_lsame(type, 'H')) {
		itype = 3;
	} else if (mc_blas_lsame(type, 'B')) {
		itype = 4;
	} else if (mc_blas_lsame(type, 'Q')) {
		itype = 5;
	} else if (mc_blas_lsame(type, 'Z')) {
		itype = 6;
	} else {
		itype = -1;
	}

	if (itype == -1) {
		*info = -1;
	} else if (cfrom == zero || mc_isnan(cfrom)) {
		*info = -4;
	} else if (mc_isnan(cto)) {
		*info = -5;
	} else if (m < 0) {
		*info = -6;
	} else if (n < 0 || (itype == 4 && n != m) || (itype == 5 && n != m)) {
		*info = -7;
	} else if (itype <= 3 && lda < mc_maxmag(1, m)) {
		*info = -9;
	} else if (itype >= 4) {
		if (kl < 0 || kl > mc_maxmag(m - 1, 0)) {
			*info = -2;
		} else {
			if (
				   (ku < 0 || ku > mc_maxmag(n - 1, 0))
				|| ((itype == 4 || itype == 5) && kl != ku)
			) {
				*info = -3;
			} else if (
				   (itype == 4 && lda < (kl + 1))
				|| (itype == 5 && lda < (ku + 1))
				|| (itype == 6 && lda < (2 * kl + ku + 1))
			) {
				*info = -9;
			}
		}
	}

	if (*info != 0) {
		mc_blas_xerbla("DLASCL", -(*info));
		return;
	}

	if (n == 0 || m == 0) {
		return;
	}

	smlnum = mc_lapack_dlamch('S');
	bignum = one / smlnum;

	cfromc = cfrom;
	ctoc   = cto;

F10:
	cfrom1 = cfromc * smlnum;
	if (cfrom1 == cfromc) {
		mul  = ctoc / cfromc;
		done = 1;
		cto1 = ctoc;
		mc_unused(cto1);
	} else {
		cto1 = ctoc / bignum;
		if (cto1 == ctoc) {
			mul    = ctoc;
			done   = 1;
			cfromc = one;
		} else if (mc_fabs(cfrom1) > mc_fabs(ctoc) && ctoc != zero) {
			mul    = smlnum;
			done   = 0;
			cfromc = cfrom1;
		} else if (mc_fabs(cto1) > mc_fabs(cfromc)) {
			mul  = bignum;
			done = 0;
			ctoc = cto1;
		} else {
			mul  = ctoc / cfromc;
			done = 1;
		}
	}

	if (itype == 0) {
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= m; ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 1) {
		for (j = 1; j <= n; ++j) {
			for (i = j; i <= m; ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 2) {
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= mc_minmag(j, m); ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 3) {
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= mc_minmag(j + 1, m); ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 4) {
		k3 = kl + 1;
		k4 = n + 1;
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= mc_minmag(k3, k4 - j); ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 5) {
		k1 = ku + 2;
		k3 = ku + 1;
		for (j = 1; j <= n; ++j) {
			for (i = mc_maxmag(k1 - j, 1); i <= k3; ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 6) {
		k1 = kl + ku + 2;
		k2 = kl + 1;
		k3 = (2 * kl) + ku + 1;
		k4 = kl + ku + 1 + m;
		for (j = 1; j <= n; ++j) {
			for (i = mc_maxmag(k1 - j, k2); i <=  mc_minmag(k3, k4 - j); ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	}
	if (!done) {
		goto F10;
	}
}

#pragma mark - mc_lapack_llascl -

MC_TARGET_FUNC void mc_lapack_llascl(const char type, const int kl, const int ku, long double cfrom, long double cto, const int m, const int n, long double * a, const int lda, int * info)
{
	const long double one = 1.0L, zero = 0.0L;

	int done;
	int i, itype, j, k1, k2, k3, k4;
	long double bignum, cfrom1, cfromc, cto1, ctoc, mul, smlnum;

	*info = 0;

	if (mc_blas_lsame(type, 'G')) {
		itype = 0;
	} else if (mc_blas_lsame(type, 'L')) {
		itype = 1;
	} else if (mc_blas_lsame(type, 'U')) {
		itype = 2;
	} else if (mc_blas_lsame(type, 'H')) {
		itype = 3;
	} else if (mc_blas_lsame(type, 'B')) {
		itype = 4;
	} else if (mc_blas_lsame(type, 'Q')) {
		itype = 5;
	} else if (mc_blas_lsame(type, 'Z')) {
		itype = 6;
	} else {
		itype = -1;
	}

	if (itype == -1) {
		*info = -1;
	} else if (cfrom == zero || mc_isnan(cfrom)) {
		*info = -4;
	} else if (mc_isnan(cto)) {
		*info = -5;
	} else if (m < 0) {
		*info = -6;
	} else if (n < 0 || (itype == 4 && n != m) || (itype == 5 && n != m)) {
		*info = -7;
	} else if (itype <= 3 && lda < mc_maxmag(1, m)) {
		*info = -9;
	} else if (itype >= 4) {
		if (kl < 0 || kl > mc_maxmag(m - 1, 0)) {
			*info = -2;
		} else {
			if (
				   (ku < 0 || ku > mc_maxmag(n - 1, 0))
				|| ((itype == 4 || itype == 5) && kl != ku)
			) {
				*info = -3;
			} else if (
				   (itype == 4 && lda < (kl + 1))
				|| (itype == 5 && lda < (ku + 1))
				|| (itype == 6 && lda < (2 * kl + ku + 1))
			) {
				*info = -9;
			}
		}
	}

	if (*info != 0) {
		mc_blas_xerbla("LLASCL", -(*info));
		return;
	}

	if (n == 0 || m == 0) {
		return;
	}

	smlnum = mc_lapack_llamch('S');
	bignum = one / smlnum;

	cfromc = cfrom;
	ctoc   = cto;

F10:
	cfrom1 = cfromc * smlnum;
	if (cfrom1 == cfromc) {
		mul  = ctoc / cfromc;
		done = 1;
		cto1 = ctoc;
		mc_unused(cto1);
	} else {
		cto1 = ctoc / bignum;
		if (cto1 == ctoc) {
			mul    = ctoc;
			done   = 1;
			cfromc = one;
		} else if (mc_fabsl(cfrom1) > mc_fabsl(ctoc) && ctoc != zero) {
			mul    = smlnum;
			done   = 0;
			cfromc = cfrom1;
		} else if (mc_fabsl(cto1) > mc_fabsl(cfromc)) {
			mul  = bignum;
			done = 0;
			ctoc = cto1;
		} else {
			mul  = ctoc / cfromc;
			done = 1;
		}
	}

	if (itype == 0) {
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= m; ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 1) {
		for (j = 1; j <= n; ++j) {
			for (i = j; i <= m; ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 2) {
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= mc_minmag(j, m); ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 3) {
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= mc_minmag(j + 1, m); ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 4) {
		k3 = kl + 1;
		k4 = n + 1;
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= mc_minmag(k3, k4 - j); ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 5) {
		k1 = ku + 2;
		k3 = ku + 1;
		for (j = 1; j <= n; ++j) {
			for (i = mc_maxmag(k1 - j, 1); i <= k3; ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	} else if (itype == 6) {
		k1 = kl + ku + 2;
		k2 = kl + 1;
		k3 = (2 * kl) + ku + 1;
		k4 = kl + ku + 1 + m;
		for (j = 1; j <= n; ++j) {
			for (i = mc_maxmag(k1 - j, k2); i <=  mc_minmag(k3, k4 - j); ++i) {
				mc_blas_matrix_at(a, lda, n, i, j) = mc_blas_matrix_at(a, lda, n, i, j) * mul;
			}
		}
	}
	if (!done) {
		goto F10;
	}
}

#endif /* !MC_LAPACKE_LASCL_H */

/* EOF */
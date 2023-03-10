//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lsamen.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>

#ifndef MC_LAPACKE_LSAMEN_H
#define MC_LAPACKE_LSAMEN_H

#pragma mark - mc_lapack_lsamen -

MC_TARGET_FUNC int mc_lapack_lsamen(const int n, const char * ca, const char * cb)
{
	int i, lsamen;
	lsamen = 0;
#	if MC_TARGET_CPP98
	if (!(n < 1 || ::strlen(ca) < mc_cast(size_t, n) || ::strlen(cb) < mc_cast(size_t, n))) {
		lsamen = 1;
		for (i = 0; i < n; ++i) {
			if (!mc_blas_lsame(ca[i], cb[i])) {
				lsamen = 0;
				break;
			}
		}
	}
#	else
	if (!(n < 1 || strlen(ca) < mc_cast(size_t, n) || strlen(cb) < mc_cast(size_t, n))) {
		lsamen = 1;
		for (i = 0; i < n; ++i) {
			if (!mc_blas_lsame(ca[i], cb[i])) {
				lsamen = 0;
				break;
			}
		}
	}
#	endif
	return lsamen;
}

#endif /* !MC_LAPACKE_LSAMEN_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_ilauplo.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>

#ifndef MC_LAPACKE_ILAUPLO_H
#define MC_LAPACKE_ILAUPLO_H

#pragma mark - mc_lapack_ilauplo -

MC_TARGET_FUNC int mc_lapack_ilauplo(const char uplo)
{
	int ilauplo;

	switch (uplo)
	{
		case 'U':
		case 'u':
			ilauplo = 121;
		break;
		case 'L':
		case 'l':
			ilauplo = 122;
		break;
		default:
			ilauplo = -1;
	}
	return ilauplo;
}

#endif /* !MC_LAPACKE_ILAUPLO_H */

/* EOF */
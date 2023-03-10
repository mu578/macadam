//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_iladiag.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>

#ifndef MC_LAPACKE_ILADIAG_H
#define MC_LAPACKE_ILADIAG_H

#pragma mark - mc_lapack_iladiag -

MC_TARGET_FUNC int mc_lapack_iladiag(const char diag)
{
	int iladiag;

	switch (diag)
	{
		case 'N':
		case 'n':
			iladiag = 131;
		break;
		case 'U':
		case 'u':
			iladiag = 132;
		break;
		default:
			iladiag = -1;
	}
	return iladiag;
}

#endif /* !MC_LAPACKE_ILADIAG_H */

/* EOF */
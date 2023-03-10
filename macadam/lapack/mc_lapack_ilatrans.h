//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_ilatrans.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>

#ifndef MC_LAPACKE_ILATRANS_H
#define MC_LAPACKE_ILATRANS_H

#pragma mark - mc_lapack_ilatrans -

MC_TARGET_FUNC int mc_lapack_ilatrans(const char trans)
{
	int ilatrans;

	switch (trans)
	{
		case 'N':
		case 'n':
			ilatrans = 111;
		break;
		case 'T':
		case 't':
			ilatrans = 112;
		break;
		case 'C':
		case 'c':
			ilatrans = 113;
		break;
		default:
			ilatrans = -1;
	}
	return ilatrans;
}

#endif /* !MC_LAPACKE_ILATRANS_H */

/* EOF */
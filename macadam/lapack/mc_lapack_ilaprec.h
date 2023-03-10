//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_ilaprec.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>

#ifndef MC_LAPACKE_ILAPREC_H
#define MC_LAPACKE_ILAPREC_H

#pragma mark - mc_lapack_ilaprec -

MC_TARGET_FUNC int mc_lapack_ilaprec(const char prec)
{
	int ilaprec;

	switch (prec)
	{
		case 'S':
		case 's':
			ilaprec = 211;
		break;
		case 'D':
		case 'd':
			ilaprec = 212;
		break;
		case 'I':
		case 'i':
			ilaprec = 213;
		break;
		case 'X':
		case 'x':
		case 'E':
		case 'e':
			ilaprec = 214;
		break;
		default:
			ilaprec = -1;
	}
	return ilaprec;
}

#endif /* !MC_LAPACKE_ILAPREC_H */

/* EOF */
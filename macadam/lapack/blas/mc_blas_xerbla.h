//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_blas_xerbla.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_BLAS_XERBLA_H
#define MC_BLAS_XERBLA_H

#pragma mark - mc_blas_xerbla -

MC_TARGET_FUNC void mc_blas_xerbla(const char * srname, const int info)
{
#	if MC_TARGET_CPP98
	::fprintf(stderr,
		"** On entry to %8s, parameter number %2d had an illegal value.\n"
		, srname
		, info
	);
#	else
	fprintf(stderr,
		"** On entry to %8s, parameter number %2d had an illegal value.\n"
		, srname
		, info
	);
#	endif
}

#endif /* !MC_BLAS_XERBLA_H */

/* EOF */
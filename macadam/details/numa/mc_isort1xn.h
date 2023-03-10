//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_isort1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_icopy1xn.h>
#include <macadam/details/numa/mc_sort1xn.h>

#ifndef MC_ISORT1XN_H 
#define MC_ISORT1XN_H 

#pragma mark - mc_isort1xn -

MC_TARGET_FUNC void mc_isort1xn(const int n, int * y, const int * x, const int f)
{
	if (x != y) {
		mc_icopy1xn(n, y, x);
	}
	mc_sort1xn_type(int, n, y, f);
}

#endif /* !MC_ISORT1XN_H  */

/* EOF */
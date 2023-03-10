//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_icopy1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/details/mc_mem.h>

#ifndef MC_ICOPY1XN_H
#define MC_ICOPY1XN_H

#pragma mark - mc_icopy1xn -

MC_TARGET_FUNC void mc_icopy1xn(const int n, int * y, const int * x)
{
	mc_base_memcpy_type(int, n, y, x);
}

#endif /* !MC_ICOPY1XN_H */

/* EOF */
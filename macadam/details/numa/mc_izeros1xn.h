//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_izeros1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/details/mc_mem.h>

#ifndef MC_IZEROS1XN_H
#define MC_IZEROS1XN_H

#pragma mark - mc_izeros1xn -

MC_TARGET_FUNC void mc_izeros1xn(const int n, int * x)
{
	mc_base_memzero_type(int, n, x);
}

#endif /* !MC_IZEROS1XN_H */

/* EOF */
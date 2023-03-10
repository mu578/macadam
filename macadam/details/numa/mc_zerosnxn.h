//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zerosnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_zeros1xn.h>

#ifndef MC_ZEROSNXN_H
#define MC_ZEROSNXN_H

#pragma mark - mc_zerosnxn -

MC_TARGET_FUNC void mc_zerosnxnf(const int n, float * a)
{
	mc_zeros1xnf(n * n, a);
}

MC_TARGET_FUNC void mc_zerosnxn(const int n, double * a)
{
	mc_zeros1xn(n * n, a);
}

MC_TARGET_FUNC void mc_zerosnxnl(const int n, long  double * a)
{
	mc_zeros1xnl(n * n, a);
}

#endif /* !MC_ZEROSNXN_H */

/* EOF */
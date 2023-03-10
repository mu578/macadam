//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zerosmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_zeros1xn.h>

#ifndef MC_ZEROSMXN_H
#define MC_ZEROSMXN_H

#pragma mark - mc_zerosmxn -

MC_TARGET_FUNC void mc_zerosmxnf(const int m, const int n, float * a)
{
	mc_zeros1xnf(m * n, a);
}

MC_TARGET_FUNC void mc_zerosmxn(const int m, const int n, double * a)
{
	mc_zeros1xn(m * n, a);
}

MC_TARGET_FUNC void mc_zerosmxnl(const int m, const int n, long  double * a)
{
	mc_zeros1xnl(m * n, a);
}

#endif /* !MC_ZEROSMXN_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_onesnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_ones1xn.h>

#ifndef MC_ONESNXN_H
#define MC_ONESNXN_H

#pragma mark - mc_onesnxn -

MC_TARGET_FUNC void mc_onesnxnf(const int n, float * a)
{
	mc_ones1xnf(n * n, a);
}

MC_TARGET_FUNC void mc_onesnxn(const int n, double * a)
{
	mc_ones1xn(n * n, a);
}

MC_TARGET_FUNC void mc_onesnxnl(const int n, long  double * a)
{
	mc_ones1xnl(n * n, a);
}

#endif /* !MC_ONESNXN_H */

/* EOF */
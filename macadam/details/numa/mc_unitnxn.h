//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_unitnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_unitmxn.h>

#ifndef MC_UNITNXN_H
#define MC_UNITNXN_H

#pragma mark - mc_unitnxn -

MC_TARGET_FUNC void mc_unitnxnf(const int n, float * a)
{
	mc_unitmxnf(n, n, a);
}

MC_TARGET_FUNC void mc_unitnxn(const int n, double * a)
{
	mc_unitmxn(n, n, a);
}

MC_TARGET_FUNC void mc_unitnxnl(const int n, long double * a)
{
	mc_unitmxnl(n, n, a);
}

#endif /* !MC_UNITNXN_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_onesmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_ones1xn.h>

#ifndef MC_ONESMXN_H
#define MC_ONESMXN_H

#pragma mark - mc_onesmxn -

MC_TARGET_FUNC void mc_onesmxnf(const int m, const int n, float * a)
{
	mc_ones1xnf(m * n, a);
}

MC_TARGET_FUNC void mc_onesmxn(const int m, const int n, double * a)
{
	mc_ones1xn(m * n, a);
}

MC_TARGET_FUNC void mc_onesmxnl(const int m, const int n, long  double * a)
{
	mc_ones1xnl(m * n, a);
}

#endif /* !MC_ONESMXN_H */

/* EOF */
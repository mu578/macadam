//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_copynxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_copymxn.h>

#ifndef MC_COPYNXN_H
#define MC_COPYNXN_H

#pragma mark - mc_copynxn -

MC_TARGET_FUNC void mc_copynxnf(const int n, float * b, const float * a)
{
	mc_copymxnf(n, n, b, a);
}

MC_TARGET_FUNC void mc_copynxnff(const int n, double * b, const float * a)
{
	mc_copymxnff(n, n, b, a);
}

MC_TARGET_FUNC void mc_copynxn(const int n, double * b, const double * a)
{
	mc_copymxn(n, n, b, a);
}

MC_TARGET_FUNC void mc_copynxnl(const int n, long double * b, const long double * a)
{
	mc_copymxnl(n, n, b, a);
}

#endif /* !MC_COPYNXN_H */

/* EOF */
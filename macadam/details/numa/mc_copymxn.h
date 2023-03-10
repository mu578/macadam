//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_copymxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_copy1xn.h>

#ifndef MC_COPYMXN_H
#define MC_COPYMXN_H

#pragma mark - mc_copymxn -

MC_TARGET_FUNC void mc_copymxnf(const int m, const int n, float * b, const float * a)
{
	mc_copy1xnf(m * n, b, a);
}

MC_TARGET_FUNC void mc_copymxnff(const int m, const int n, double * b, const float * a)
{
	mc_copy1xnff(m * n, b, a);
}

MC_TARGET_FUNC void mc_copymxn(const int m, const int n, double * b, const double * a)
{
	mc_copy1xn(m * n, b, a);
}

MC_TARGET_FUNC void mc_copymxnl(const int m, const int n, long double * b, const long double * a)
{
	mc_copy1xnl(m * n, b, a);
}

#endif /* !MC_COPYMXN_H */

/* EOF */
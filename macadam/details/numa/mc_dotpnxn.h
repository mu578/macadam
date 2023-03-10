//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_dotpnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_dotpmxn.h>

#ifndef MC_DOTPNXN_H
#define MC_DOTPNXN_H

#pragma mark - mc_dotpnxn -

MC_TARGET_FUNC void mc_dotpnxnf(const int n, float * MC_TARGET_RESTRICT c, const float * a, const float * b, const int d, const int f)
{
	mc_dotpmxnf(n, n, n, c, a, b, d, f);
}

MC_TARGET_FUNC void mc_dotpnxn(const int n, double * MC_TARGET_RESTRICT c, const double * a, const double * b, const int d, const int f)
{
	mc_dotpmxn(n, n, n, c, a, b, d, f);
}

MC_TARGET_FUNC void mc_dotpnxnl(const int n, long double * MC_TARGET_RESTRICT c, const long double * a, const long double * b, const int d, const int f)
{
	mc_dotpmxnl(n, n, n, c, a, b, d, f);
}

#endif /* !MC_DOTPNXN_H */

/* EOF */
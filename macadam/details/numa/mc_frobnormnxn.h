//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_frobnormnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_frobnormmxn.h>

#ifndef MC_FROBNORMNXN_H
#define MC_FROBNORMNXN_H

#pragma mark - mc_frobnormnxn -

MC_TARGET_FUNC float mc_frobnormnxnf(const int n, const float * a)
{
	return mc_frobnormmxnf(n, n, a);
}

MC_TARGET_FUNC double mc_frobnormnxnff(const int n, const float * a)
{
	return mc_frobnormmxnff(n, n, a);
}

MC_TARGET_FUNC double mc_frobnormnxn(const int n, const double * a)
{
	return mc_frobnormmxn(n, n, a);
}

MC_TARGET_FUNC long double mc_frobnormnxnl(const int n, const long double * a)
{
	return mc_frobnormmxnl(n, n, a);
}

#endif /* !MC_FROBNORMNXN_H */

/* EOF */
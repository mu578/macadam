//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_trspnxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_trspmxn.h>

#ifndef MC_TRSPNXN_H
#define MC_TRSPNXN_H

#pragma mark - mc_trspnxn -

MC_TARGET_FUNC void mc_trspnxnf(const int n, float * at, const float * a)
{
//!# Returning transpose of A.
	mc_trspmxnf(n, n, at, a);
}

MC_TARGET_FUNC void mc_trspnxn(const int n, double * at, const double * a)
{
//!# Returning transpose of A.
	mc_trspmxn(n, n, at, a);
}

MC_TARGET_FUNC void mc_trspnxnl(const int n, long double * at, const long double * a)
{
//!# Returning transpose of A.
	mc_trspmxnl(n, n, at, a);
}

#endif /* !MC_TRSPNXN_H */

/* EOF */
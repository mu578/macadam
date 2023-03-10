//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_frobnormmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_l2norm1xn.h>

#ifndef MC_FROBNORMMXN_H
#define MC_FROBNORMMXN_H

#pragma mark - mc_frobnormmxn -

MC_TARGET_FUNC float mc_frobnormmxnf(const int m, const int n, const float * a)
{
	return mc_l2norm1xnf(m * n, a);
}

MC_TARGET_FUNC double mc_frobnormmxnff(const int m, const int n, const float * a)
{
	return mc_l2norm1xnff(m * n, a);
}

MC_TARGET_FUNC double mc_frobnormmxn(const int m, const int n, const double * a)
{
	return mc_l2norm1xn(m * n, a);
}

MC_TARGET_FUNC long double mc_frobnormmxnl(const int m, const int n, const long double * a)
{
	return mc_l2norm1xnl(m * n, a);
}

#endif /* !MC_FROBNORMMXN_H */

/* EOF */
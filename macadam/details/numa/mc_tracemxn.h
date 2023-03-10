//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_tracemxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_minmag.h>
#include <macadam/details/numa/mc_tracenxn.h>

#ifndef MC_TRACEMXN_H
#define MC_TRACEMXN_H

#pragma mark - mc_tracemxn -

MC_TARGET_FUNC float mc_tracemxnf(const int m, const int n, const float * a)
{
	if (m >= n) {
		return mc_tracenxnf(n, a);
	}
	return 0.0f;
}

MC_TARGET_FUNC double mc_tracemxnff(const int m, const int n, const float * a)
{
	if (m >= n) {
		return mc_tracenxnff(n, a);
	}
	return 0.0;
}

MC_TARGET_FUNC double mc_tracemxn(const int m, const int n, const double * a)
{
	if (m >= n) {
		return mc_tracenxn(n, a);
	}
	return 0.0;
}

MC_TARGET_FUNC long double mc_tracemxnl(const int m, const int n, const long double * a)
{
	if (m >= n) {
		return mc_tracenxnl(n, a);
	}
	return 0.0L;
}

#endif /* !MC_TRACEMXN_H */

/* EOF */
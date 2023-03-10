//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ones1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_ONES1XN_H
#define MC_ONES1XN_H

#pragma mark - mc_ones1xn -

MC_TARGET_FUNC void mc_ones1xnf(const int n, float * x)
{
#	if MC_TARGET_CPP98
	::std::fill_n(x, n, 1.0f);
#	else
	int i = 0;
	for (; i < n; i++) {
		x[i] = 1.0f;
	}
#	endif
}

MC_TARGET_FUNC void mc_ones1xn(const int n, double * x)
{
#	if MC_TARGET_CPP98
	::std::fill_n(x, n, 1.0);
#	else
	int i = 0;
	for (; i < n; i++) {
		x[i] = 1.0;
	}
#	endif
}

MC_TARGET_FUNC void mc_ones1xnl(const int n, long double * x)
{
#	if MC_TARGET_CPP98
	::std::fill_n(x, n, 1.0L);
#	else
	int i = 0;
	for (; i < n; i++) {
		x[i] = 1.0L;
	}
#	endif
}

#endif /* !MC_ONES1XN_H */

/* EOF */
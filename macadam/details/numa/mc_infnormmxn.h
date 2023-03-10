//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_infnormmxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_a2sum1xn.h>

#ifndef MC_INFNORMMXN_H
#define MC_INFNORMMXN_H

#pragma mark - mc_infnormmxn -

MC_TARGET_FUNC float mc_infnormmxnf(const int m, const int n, const float * a, const int f)
{
//!# Requires a[m x n]. Returning the infinity norm of a.
//!# f=0: computing the maximum of the absolute row sums.
//!# f=1: computing the minimum of the absolute row sums.
	int i = 0;
	float nrm = 0.0f, sum;
	for (; i < m ; i++) {
		const float * x = a + (n * i);
		sum             = mc_a2sum1xnf(n, x);
		nrm             = (f == 1) ? ((nrm > sum) ? sum : nrm) : ((nrm > sum) ? nrm : sum);
	}
	return nrm;
}

MC_TARGET_FUNC double mc_infnormmxnff(const int m, const int n, const float * a, const int f)
{
//!# Requires a[m x n]. Returning the infinity norm of a.
//!# f=0: computing the maximum of the absolute row sums.
//!# f=1: computing the minimum of the absolute row sums.
	int i = 0;
	double nrm = 0.0, sum;
	for (; i < m ; i++) {
		const float * x = a + (n * i);
		sum             = mc_a2sum1xnff(n, x);
		nrm             = (f == 1) ? ((nrm > sum) ? sum : nrm) : ((nrm > sum) ? nrm : sum);
	}
	return nrm;
}

MC_TARGET_FUNC double mc_infnormmxn(const int m, const int n, const double * a, const int f)
{
//!# Requires a[m x n]. Returning the infinity norm of a.
//!# f=0: computing the maximum of the absolute row sums.
//!# f=1: computing the minimum of the absolute row sums.
	int i = 0;
	double nrm = 0.0, sum;
	for (; i < m ; i++) {
		const double * x = a + (n * i);
		sum              = mc_a2sum1xn(n, x);
		nrm              = (f == 1) ? ((nrm > sum) ? sum : nrm) : ((nrm > sum) ? nrm : sum);
	}
	return nrm;
}

MC_TARGET_FUNC long double mc_infnormmxnl(const int m, const int n, const long double * a, const int f)
{
//!# Requires a[m x n]. Returning the infinity norm of a.
//!# f=0: computing the maximum of the absolute row sums.
//!# f=1: computing the minimum of the absolute row sums.
	int i = 0;
	long double nrm = 0.0L, sum;
	for (; i < m ; i++) {
		const long double * x = a + (n * i);
		sum                   = mc_a2sum1xnl(n, x);
		nrm                   = (f == 1) ? ((nrm > sum) ? sum : nrm) : ((nrm > sum) ? nrm : sum);
	}
	return nrm;
}

#endif /* !MC_INFNORMMXN_H */

/* EOF */
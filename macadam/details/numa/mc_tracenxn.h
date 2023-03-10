//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_tracenxn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_TRACENXN_H
#define MC_TRACENXN_H

#pragma mark - mc_tracenxn -

MC_TARGET_FUNC float mc_tracenxnf(const int n, const float * a)
{
	int i      = 0;
	float diag = 0.0f;
	for (; i < n; i++) {
		diag = diag + a[(n * i) + i];
	}
	return diag;
}

MC_TARGET_FUNC double mc_tracenxnff(const int n, const float * a)
{
	int i       = 0;
	double diag = 0.0;
	for (; i < n; i++) {
		diag = diag + mc_cast(double, a[(n * i) + i]);
	}
	return diag;
}

MC_TARGET_FUNC double mc_tracenxn(const int n, const double * a)
{
	int i       = 0;
	double diag = 0.0;
	for (; i < n; i++) {
		diag = diag + a[(n * i) + i];
	}
	return diag;
}

MC_TARGET_FUNC long double mc_tracenxnl(const int n, const long double * a)
{
	int i            = 0;
	long double diag = 0.0;
	for (; i < n; i++) {
		diag = diag + a[(n * i) + i];
	}
	return diag;
}

#endif /* !MC_TRACENXN_H */

/* EOF */
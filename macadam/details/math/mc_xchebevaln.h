//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_xchebevaln.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_XCHEBEVALN_H
#define MC_XCHEBEVALN_H

#pragma mark - mc_xchebevaln -

MC_TARGET_PROC float mc_xchebevalnf(float x, const float * a, const unsigned int n)
{
//!# Evaluating Chebyshev sum using Clenshaw algorithm.
	int i, m;
	float s = 0.0f, b0 = 0.0f, b1 = 0.0f, b2 = 0.0f;
	if (mc_nonnullptr(a) && (n > 0 && n < 256)) {
		m = mc_cast(int, n);
		x = 2.0f * x;
		for (i = m - 1; i >= 0; i--) {
			b2 = b1;
			b1 = b0;
			b0 = x * b1 - b2 + a[i];
		}
		s = 0.5f * (b0 - b2);
	}
	return s;
}

MC_TARGET_PROC double mc_xchebevaln(double x, const double * a, const unsigned int n)
{
//!# Evaluating Chebyshev sum using Clenshaw algorithm.
	int i, m;
	double s = 0.0, b0 = 0.0, b1 = 0.0, b2 = 0.0;
	if (mc_nonnullptr(a) && (n > 0 && n < 256)) {
		m = mc_cast(int, n);
		x = 2.0 * x;
		for (i = m - 1; i >= 0; i--) {
			b2 = b1;
			b1 = b0;
			b0 = x * b1 - b2 + a[i];
		}
		s = 0.5 * (b0 - b2);
	}
	return s;
}

MC_TARGET_PROC long double mc_xchebevalnl(long double x, const long double * a, const unsigned int n)
{
//!# Evaluating Chebyshev sum using Clenshaw algorithm.
	int i, m;
	long double s = 0.0L, b0 = 0.0L, b1 = 0.0L, b2 = 0.0L;
	if (mc_nonnullptr(a) && (n > 0 && n < 256)) {
		m = mc_cast(int, n);
		x = 2.0L * x;
		for (i = m - 1; i >= 0; i--) {
			b2 = b1;
			b1 = b0;
			b0 = x * b1 - b2 + a[i];
		}
		s = 0.5L * (b0 - b2);
	}
	return s;
}

#endif /* !MC_XCHEBEVALN_H */

/* EOF */
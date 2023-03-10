//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_chbevl.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_CHBEVL_H
#define MC_CHBEVL_H

#pragma mark - mc_chbevl -

MC_TARGET_PROC float mc_chbevlf(const float x, const float * a, const unsigned int n)
{
	int i;
	float s = 0.0f, b0 = 0.0f, b1 = 0.0f, b2 = 0.0f;
	const float * p;
	if (mc_nonnullptr(a) && (n > 0 && n < 256)) {
		p  = a;
		b0 = *p++;
		i  = mc_cast(int, n) - 1;
		do {
			b2 = b1;
			b1 = b0;
			b0 = x * b1 -  b2 + *p++;
		} while (--i);
		s = 0.5f * (b0 - b2);
	}
	return s;
}

MC_TARGET_PROC double mc_chbevl(const double x, const double * a, const unsigned int n)
{
	int i;
	double s = 0.0, b0 = 0.0, b1 = 0.0, b2 = 0.0;
	const double * p;
	if (mc_nonnullptr(a) && (n > 0 && n < 256)) {
		p  = a;
		b0 = *p++;
		i  = mc_cast(int, n) - 1;
		do {
			b2 = b1;
			b1 = b0;
			b0 = x * b1 -  b2 + *p++;
		} while (--i);
		s = 0.5 * (b0 - b2);
	}
	return s;
}

MC_TARGET_PROC long double mc_chbevll(const long double x, const long double * a, const unsigned int n)
{
	int i;
	long double s = 0.0L, b0 = 0.0L, b1 = 0.0L, b2 = 0.0L;
	const long double * p;
	if (mc_nonnullptr(a) && (n > 0 && n < 256)) {
		p  = a;
		b0 = *p++;
		i  = mc_cast(int, n) - 1;
		do {
			b2 = b1;
			b1 = b0;
			b0 = x * b1 -  b2 + *p++;
		} while (--i);
		s = 0.5L * (b0 - b2);
	}
	return s;
}

#endif /* !MC_CHBEVL_H */

/* EOF */
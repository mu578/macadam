//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_hermite_hen.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>

#ifndef MC_HERMITE_HEN_H
#define MC_HERMITE_HEN_H

#pragma mark - mc_hermite_hen -

MC_TARGET_FUNC float mc_hermite_henf(const unsigned int n, const float x)
{
//!# NTH probabilists' Hermite polynomial.
	float hi       = 0.0f;
	float h0       = 0.0f;
	float h1       = 0.0f;
	unsigned int i = 1;
	if (mc_isnan(x) || mc_isinf(x)) {
		return MCK_NAN;
	}
	if (n >= 1 && n < MCLIMITS_IMAX) {
		h0 = 1.0f;
		h1 = x;
		for (; i < n; i++) {
			hi = x * h1 - mc_cast(float, i) * h0;
			h0 = h1;
			h1 = hi;
		}
	}
	return h1;
}

MC_TARGET_FUNC double mc_hermite_hen(const unsigned int n, const double x)
{
//!# NTH probabilists' Hermite polynomial.
	double hi      = 0.0;
	double h0      = 0.0;
	double h1      = 0.0;
	unsigned int i = 1;
	if (mc_isnan(x) || mc_isinf(x)) {
		return MCK_NAN;
	}
	if (n >= 1 && n < MCLIMITS_IMAX) {
		h0 = 1.0;
		h1 = x;
		for (; i < n; i++) {
			hi = x * h1 - mc_cast(double, i) * h0;
			h0 = h1;
			h1 = hi;
		}
	}
	return h1;
}

MC_TARGET_FUNC long double mc_hermite_henl(const unsigned int n, const long double x)
{
//!# NTH probabilists' Hermite polynomial.
	long double hi = 0.0L;
	long double h0 = 0.0L;
	long double h1 = 0.0L;
	unsigned int i = 1;
	if (mc_isnan(x) || mc_isinf(x)) {
		return MCK_NAN;
	}
	if (n >= 1 && n < MCLIMITS_IMAX) {
		h0 = 1.0L;
		h1 = x;
		for (; i < n; i++) {
			hi = x * h1 - mc_cast(long double, i) * h0;
			h0 = h1;
			h1 = hi;
		}
	}
	return h1;
}

#endif /* !MC_HERMITE_HEN_H */

/* EOF */
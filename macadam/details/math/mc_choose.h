//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_choose.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_beta.h>
#include <macadam/details/math/mc_ceil.h>
#include <macadam/details/math/mc_factorial.h>

#ifndef MC_CHOOSE_H
#define MC_CHOOSE_H

#pragma mark - mc_choose -

MC_TARGET_FUNC float mc_choosef(const unsigned int n, const unsigned int k)
{
	float r                          = MCLIMITS_MAXF;
	const unsigned int factorial_max = 35U;

	if (n < MCLIMITS_UIMAX && k < MCLIMITS_UIMAX) {
		if (k > n) {
			return MCLIMITS_MAXF;
		} else if ((k == 0) || (k == n)) {
			return 1.0f;
		} else if ((k == 1) || (k == n - 1)) {
			return n < 0x1000001U ? mc_cast(float, n) : MCLIMITS_MAXF;
		}
		if (n < factorial_max) {
			r = mc_factorialf(n);
			r = r / mc_factorialf(n - k);
			r = r / mc_factorialf(k);
			r = mc_ceilf(r - 0.5f);
		} else {
			if (n < 0x1000001U && k < 0x1000001U) {
				r = (k < (n - k))
					? mc_cast(float, (k - 0)) * mc_betaf(mc_cast(float, (k + 0)), mc_cast(float, (n - k + 1)))
					: mc_cast(float, (n - k)) * mc_betaf(mc_cast(float, (k + 1)), mc_cast(float, (n - k + 0)))
				;
			} else {
				r = mc_cast_expr(float, (k < (n - k))
					? mc_cast(double, (k - 0)) * mc_beta(mc_cast(double, (k + 0)), mc_cast(double, (n - k + 1)))
					: mc_cast(double, (n - k)) * mc_beta(mc_cast(double, (k + 1)), mc_cast(double, (n - k + 0)))
				);
			}
			if (r != 0) {
				r = 1.0f / r;
				r = mc_ceilf(r - 0.5f);
			} else {
				r = MCLIMITS_MAXF;
			}
		}
	}
	return r;
}

MC_TARGET_FUNC double mc_choose(const unsigned int n, const unsigned int k)
{
	double r                         = MCLIMITS_MAX;
	const unsigned int factorial_max = 171U;
	if (n < MCLIMITS_UIMAX && k < MCLIMITS_UIMAX) {
		if (k > n) {
			return MCLIMITS_MAX;
		} else if ((k == 0) || (k == n)) {
			return 1.0;
		} else if ((k == 1) || (k == n - 1)) {
			return mc_cast(double, n);
		}
		if (n < factorial_max) {
			r = mc_factorial(n);
			r = r / mc_factorial(n - k);
			r = r / mc_factorial(k);
			r = mc_ceil(r - 0.5);
		} else {
			r = (k < (n - k))
				? mc_cast(double, (k - 0)) * mc_beta(mc_cast(double, (k + 0)), mc_cast(double, (n - k + 1)))
				: mc_cast(double, (n - k)) * mc_beta(mc_cast(double, (k + 1)), mc_cast(double, (n - k + 0)))
			;
			if (r != 0) {
				r = 1.0f / r;
				r = mc_ceil(r - 0.5);
			} else {
				r = MCLIMITS_MAX;
			}
		}
	}
	return r;
}

MC_TARGET_FUNC long double mc_choosel(const unsigned int n, const unsigned int k)
{
	long double r                    = MCLIMITS_MAXL;
#	if MC_TARGET_HAVE_LONG_DOUBLE
	const unsigned int factorial_max = 1755U;
#	else
	const unsigned int factorial_max = 171U;
#	endif
	if (n < MCLIMITS_UIMAX && k < MCLIMITS_UIMAX) {
		if (k > n) {
			return MCLIMITS_MAXL;
		} else if ((k == 0) || (k == n)) {
			return 1.0L;
		} else if ((k == 1) || (k == n - 1)) {
			return mc_cast(long double, n);
		}
		if (n < factorial_max) {
			r = mc_cast(long double, mc_factoriall(n));
			r = r / mc_cast(const long double, mc_factoriall(n - k));
			r = r / mc_cast(const long double, mc_factoriall(k));
			r = mc_ceill(r - 0.5L);
		} else {
			r = (k < (n - k))
				? mc_cast(long double, (k - 0)) * mc_betal(mc_cast(long double, (k + 0)), mc_cast(long double, (n - k + 1)))
				: mc_cast(long double, (n - k)) * mc_betal(mc_cast(long double, (k + 1)), mc_cast(long double, (n - k + 0)))
			;
			if (r != 0) {
				r = 1.0f / r;
				r = mc_ceill(r - 0.5L);
			} else {
				r = MCLIMITS_MAXL;
			}
		}
	}
	return r;
}

#endif /* !MC_CHOOSE_H */

/* EOF */
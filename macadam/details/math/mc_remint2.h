//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_remint2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fisint.h>
#include <macadam/details/math/mc_fisodd.h>
#include <macadam/details/math/mc_fisval.h>
#include <macadam/details/math/mc_floor.h>
#include <macadam/details/math/mc_itrunc64.h>

#ifndef MC_REMINT2_H
#define MC_REMINT2_H

#pragma mark - mc_remint2 -

MC_TARGET_PROC int64_t mc_remint2f(float x, float * z)
{
//!# Reduction of x=z*PI=x*PI - n*PI/2 where |z|<=1/4, result = n mod 8.
	int64_t r         = 0;
	float w;
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	const int64_t max = INT64_MAX;
#	else
	const int64_t max = mc_cast(const int64_t, 0x7FFFFFFFFFFFFFFF);
#	endif
	*z                = 0.0f;
	if (mc_fisvalf(x)) {
		if (mc_fabsf(x) <= 0.25f) {
			 r = 0;
			*z = x;
		} else if (mc_fisintf(x)) {
			if (mc_fisoddf(x, 0)) {
				r = 2;
			}
		} else {
			x = 0.5f * x;
			x = 2.0f * (x - mc_floorf(x));
			w = mc_floorf(4.0f * x);
			r = mc_itrunc64f(w - 16.0f * mc_floorf(w / 16.0f));
			if (r < (max - 1)) {
				if (r & 1) {
					++r;
					w = w + 1.0f;
				}
				 r = (mc_cast(uint64_t, r) >> 1) & 7;
				*z = x - 0.25f * w;
			} else {
				 r = 0;
				*z = x;
			}
		}
	}
	return r;
}

MC_TARGET_PROC int64_t mc_remint2(double x, double * z)
{
//!# Reduction of x=z*PI=x*PI - n*PI/2 where |z|<=1/4, result = n mod 8.
	int64_t r         = 0;
	double w;
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	const int64_t max = INT64_MAX;
#	else
	const int64_t max = mc_cast(const int64_t, 0x7FFFFFFFFFFFFFFF);
#	endif
	*z                = 0.0;
	if (mc_fisval(x)) {
		if (mc_fabs(x) <= 0.25) {
			 r = 0;
			*z = x;
		} else if (mc_fisint(x)) {
			if (mc_fisodd(x, 0)) {
				r = 2;
			}
		} else {
			x = 0.5 * x;
			x = 2.0 * (x - mc_floor(x));
			w = mc_floor(4.0 * x);
			r = mc_itrunc64(w - 16.0 * mc_floor(w / 16.0));
			if (r < (max - 1)) {
				if (r & 1) {
					++r;
					w = w + 1.0;
				}
				 r = (mc_cast(uint64_t, r) >> 1) & 7;
				*z = x - 0.25 * w;
			} else {
				 r = 0;
				*z = x;
			}
		}
	}
	return r;
}

MC_TARGET_PROC int64_t mc_remint2l(long double x, long double * z)
{
#	if MC_TARGET_HAVE_LONG_DOUBLE
//!# Reduction of x=z*PI=x*PI - n*PI/2 where |z|<=1/4, result = n mod 8.
	int64_t r         = 0;
	long double w;
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	const int64_t max = INT64_MAX;
#	else
	const int64_t max = mc_cast(const int64_t, 0x7FFFFFFFFFFFFFFF);
#	endif
	*z                = 0.0L;
	if (mc_fisvall(x)) {
		if (mc_fabsl(x) <= 0.25L) {
			 r = 0;
			*z = x;
		} else if (mc_fisintl(x)) {
			if (mc_fisoddl(x, 0)) {
				r = 2;
			}
		} else {
			x = 0.5L * x;
			x = 2.0L * (x - mc_floorl(x));
			w = mc_floorl(4.0L * x);
			r = mc_itrunc64l(w - 16.0L * mc_floorl(w / 16.0L));
			if (r < (max - 1)) {
				if (r & 1) {
					++r;
					w = w + 1.0L;
				}
				 r = (mc_cast(uint64_t, r) >> 1) & 7;
				*z = x - 0.25L * w;
			} else {
				 r = 0;
				*z = x;
			}
		}
	}
	return r;
#	else
	double y  = mc_cast(double, x), w = 0.0;
	int64_t r = mc_remint2(y, &w);
	*z        = mc_cast(long double, w);
	return r;
#	endif
}

#endif /* !MC_REMINT2_H */

/* EOF */
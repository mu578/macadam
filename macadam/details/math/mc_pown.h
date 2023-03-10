//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_pown.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_raise3.h>
#include <macadam/details/math/mc_raise4.h>
#include <macadam/details/math/mc_raise5.h>
#include <macadam/details/math/mc_raise6.h>

#ifndef MC_POWN_H
#define MC_POWN_H

#pragma mark - mc_pown -

MC_TARGET_FUNC float mc_pownf(const float x, const int y)
{
	float r;
	if (mc_isnan(x)) {
		return x;
	}
	if (mc_isinf(x)) {
		if (y < 0 && !(y & 1)) {
			return 0.0f;
		}
		if (y < 0 && y & 1) {
			return mc_copysignf(0.0f, x);
		}
		return x;
	}
	switch (y)
	{
		case 0:
			return 1.0f;
		case 1:
			return x;
		case 2:
			return mc_raise2f(x);
		case 3:
			return mc_raise3f(x);
		case 4:
			return mc_raise4f(x);
		case 5:
			return mc_raise5f(x);
		case 6:
			return mc_raise6f(x);
	}
	if (mc_iabs(y) < 0x1000001) {
		r = mc_powf(x, mc_cast(const float, y));
	} else {
		r = mc_cast(float, mc_pow(mc_cast(const double, x), mc_cast(const double, y)));
	}
	if (mc_isinf(r)) {
		r = mc_copysignf(r, x);
	}
	return r;
}

MC_TARGET_FUNC double mc_pown(const double x, const int y)
{
	double r;
	if (mc_isnan(x)) {
		return x;
	}
	if (mc_isinf(x)) {
		if (y < 0 && !(y & 1)) {
			return 0.0;
		}
		if (y < 0 && y & 1) {
			return mc_copysign(0.0, x);
		}
		return x;
	}
	switch (y)
	{
		case 0:
			return 1.0;
		case 1:
			return x;
		case 2:
			return mc_raise2(x);
		case 3:
			return mc_raise3(x);
		case 4:
			return mc_raise4(x);
		case 5:
			return mc_raise5(x);
		case 6:
			return mc_raise6(x);
	}
	r = mc_pow(x, mc_cast(const double, y));
	if (mc_isinf(r)) {
		return mc_copysign(r, x);
	}
	return r;
}

MC_TARGET_FUNC long double mc_pownl(const long double x, const int y)
{
	long double r;
	if (mc_isnan(x)) {
		return x;
	}
	if (mc_isinf(x)) {
		if (y < 0 && !(y & 1)) {
			return 0.0L;
		}
		if (y < 0 && y & 1) {
			return mc_copysignl(0.0L, x);
		}
		return x;
	}
	switch (y)
	{
		case 0:
			return 1.0L;
		case 1:
			return x;
		case 2:
			return mc_raise2l(x);
		case 3:
			return mc_raise3l(x);
		case 4:
			return mc_raise4l(x);
		case 5:
			return mc_raise5l(x);
		case 6:
			return mc_raise6l(x);
	}
	r = mc_powl(x, mc_cast(const long double, y));
	if (mc_isinf(r)) {
		return mc_copysignl(r, x);
	}
	return r;
}

#pragma mark - mc_powln -

MC_TARGET_FUNC float mc_powlnf(const float x, const long y)
{
	float r;
	if (mc_isnan(x)) {
		return x;
	}
	if (mc_isinf(x)) {
		if (y < 0 && !(y & 1)) {
			return 0.0f;
		}
		if (y < 0 && y & 1) {
			return mc_copysignf(0.0f, x);
		}
		return x;
	}
	switch (y)
	{
		case 0:
			return 1.0f;
		case 1:
			return x;
		case 2:
			return mc_raise2f(x);
		case 3:
			return mc_raise3f(x);
		case 4:
			return mc_raise4f(x);
		case 5:
			return mc_raise5f(x);
		case 6:
			return mc_raise6f(x);
	}
	if (mc_labs(y) < 0x1000001) {
		r = mc_powf(x, mc_cast(const float, y));
	} else {
#	if MC_TARGET_LONG_64BIT
		if (mc_labs(y) < 0x0020000000000000L) {
			r = mc_cast(float, mc_pow(mc_cast(const double, x), mc_cast(const double, y)));
		} else {
#		if MC_TARGET_HAVE_LONG_DOUBLE
			r = mc_cast(float, mc_powl(mc_cast(const long double, x), mc_cast(const long double, y)));
#		else
			r = MCK_INF;
#		endif
		}
#	else
		r = mc_cast(float, mc_pow(mc_cast(const double, x), mc_cast(const double, y)));
#	endif
	}
	if (mc_isinf(r)) {
		r = mc_copysignf(r, x);
	}
	return r;
}

MC_TARGET_FUNC double mc_powln(const double x, const long y)
{
	double r;
	if (mc_isnan(x)) {
		return x;
	}
	if (mc_isinf(x)) {
		if (y < 0 && !(y & 1)) {
			return 0.0;
		}
		if (y < 0 && y & 1) {
			return mc_copysign(0.0, x);
		}
		return x;
	}
	switch (y)
	{
		case 0:
			return 1.0;
		case 1:
			return x;
		case 2:
			return mc_raise2(x);
		case 3:
			return mc_raise3(x);
		case 4:
			return mc_raise4(x);
		case 5:
			return mc_raise5(x);
		case 6:
			return mc_raise6(x);
	}
#	if MC_TARGET_LONG_64BIT
	if (mc_labs(y) < 0x0020000000000000L) {
		r = mc_pow(x, mc_cast(const double, y));
	} else {
#		if MC_TARGET_HAVE_LONG_DOUBLE
		r = mc_cast(double, mc_powl(mc_cast(const long double, x), mc_cast(const long double, y)));
#		else
		r = MCK_INF;
#		endif
	}
#	else
	r = mc_pow(x, mc_cast(const double, y));
#	endif
	if (mc_isinf(r)) {
		return mc_copysign(r, x);
	}
	return r;
}

MC_TARGET_FUNC long double mc_powlnl(const long double x, const long y)
{
	long double r;
	if (mc_isnan(x)) {
		return x;
	}
	if (mc_isinf(x)) {
		if (y < 0 && !(y & 1)) {
			return 0.0L;
		}
		if (y < 0 && y & 1) {
			return mc_copysignl(0.0L, x);
		}
		return x;
	}
	switch (y)
	{
		case 0:
			return 1.0L;
		case 1:
			return x;
		case 2:
			return mc_raise2l(x);
		case 3:
			return mc_raise3l(x);
		case 4:
			return mc_raise4l(x);
		case 5:
			return mc_raise5l(x);
		case 6:
			return mc_raise6l(x);
	}
#	if !MC_TARGET_HAVE_LONG_DOUBLE && MC_TARGET_LONG_64BIT
		if (mc_labs(y) < 0x0020000000000000L) {
			r = mc_powl(x, mc_cast(const long double, y));
		} else {
			r = MCK_INF;
		}
#	else
	r = mc_powl(x, mc_cast(const long double, y));
#	endif
	if (mc_isinf(r)) {
		return mc_copysignl(r, x);
	}
	return r;
}

#pragma mark - mc_powlln -

#	if MC_TARGET_C99 || MC_TARGET_CPP11
MC_TARGET_FUNC float mc_powllnf(const float x, const long long y)
{
	float r;
	if (mc_isnan(x)) {
		return x;
	}
	if (mc_isinf(x)) {
		if (y < 0 && !(y & 1)) {
			return 0.0f;
		}
		if (y < 0 && y & 1) {
			return mc_copysignf(0.0f, x);
		}
		return x;
	}
	switch (y)
	{
		case 0:
			return 1.0f;
		case 1:
			return x;
		case 2:
			return mc_raise2f(x);
		case 3:
			return mc_raise3f(x);
		case 4:
			return mc_raise4f(x);
		case 5:
			return mc_raise5f(x);
		case 6:
			return mc_raise6f(x);
	}
	if (mc_llabs(y) < 0x1000001) {
		r = mc_powf(x, mc_cast(const float, y));
	} else if (mc_llabs(y) < 0x0020000000000000L) {
		r = mc_cast(float, mc_pow(mc_cast(const double, x), mc_cast(const double, y)));
	} else {
#	if MC_TARGET_HAVE_LONG_DOUBLE
		r = mc_cast(float, mc_powl(mc_cast(const long double, x), mc_cast(const long double, y)));
#	else
		r = MCK_INF;
#	endif
	}
	if (mc_isinf(r)) {
		r = mc_copysignf(r, x);
	}
	return r;
}

MC_TARGET_FUNC double mc_powlln(const double x, const long long y)
{
	double r;
	if (mc_isnan(x)) {
		return x;
	}
	if (mc_isinf(x)) {
		if (y < 0 && !(y & 1)) {
			return 0.0;
		}
		if (y < 0 && y & 1) {
			return mc_copysign(0.0, x);
		}
		return x;
	}
	switch (y)
	{
		case 0:
			return 1.0;
		case 1:
			return x;
		case 2:
			return mc_raise2(x);
		case 3:
			return mc_raise3(x);
		case 4:
			return mc_raise4(x);
		case 5:
			return mc_raise5(x);
		case 6:
			return mc_raise6(x);
	}
	if (mc_llabs(y) < 0x0020000000000000L) {
		r = mc_pow(x, mc_cast(double, y));
	} else {
#	if MC_TARGET_HAVE_LONG_DOUBLE
		r = mc_cast(double, mc_powl(mc_cast(long double, x), mc_cast(long double, y)));
#	else
		r = MCK_INF;
#	endif
	}
	if (mc_isinf(r)) {
		return mc_copysign(r, x);
	}
	return r;
}

MC_TARGET_FUNC long double mc_powllnl(const long double x, const long long y)
{
#	if MC_TARGET_HAVE_LONG_DOUBLE
	long double r;
	if (mc_isnan(x)) {
		return x;
	}
	if (mc_isinf(x)) {
		if (y < 0 && !(y & 1)) {
			return 0.0L;
		}
		if (y < 0 && y & 1) {
			return mc_copysignl(0.0L, x);
		}
		return x;
	}
	switch (y)
	{
		case 0:
			return 1.0L;
		case 1:
			return x;
		case 2:
			return mc_raise2l(x);
		case 3:
			return mc_raise3l(x);
		case 4:
			return mc_raise4l(x);
		case 5:
			return mc_raise5l(x);
		case 6:
			return mc_raise6l(x);
	}
	r = mc_powl(x, mc_cast(const long double, y));
	if (mc_isinf(r)) {
		return mc_copysignl(r, x);
	}
	return r;
#	else
	return mc_cast(long double, mc_powlln(mc_cast(const double, x), y));
#	endif
}
#	else
#	define mc_powllnf mc_powlnf
#	define mc_powlln  mc_powln
#	define mc_powllnl mc_polnl
#	endif

#endif /* !MC_POWN_H */

/* EOF */
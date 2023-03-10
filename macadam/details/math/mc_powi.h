//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_powi.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_raise3.h>
#include <macadam/details/math/mc_raise4.h>
#include <macadam/details/math/mc_raise5.h>
#include <macadam/details/math/mc_raise6.h>

#ifndef MC_POWI_H
#define MC_POWI_H

#pragma mark - mc_powi -

MC_TARGET_FUNC float mc_powif(const float x, const int y)
{

	if (x == 0.0f) {
		if (y == 0) {
			return 1.0f;
		} else if (y < 0) {
			return mc_copysignf(MCK_INF, x);
		}
		return ((y % 2) == 0) ? 0.0f : mc_copysignf(0.0f, x);
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
		return mc_powf(x, mc_cast(float, y));
	} else {
		return mc_cast(float, mc_pow(mc_cast(double, x), mc_cast(double, y)));
	}
}

MC_TARGET_FUNC double mc_powi(const double x, const int y)
{
	if (x == 0.0) {
		if (y == 0) {
			return 1.0;
		} else if (y < 0) {
			return mc_copysign(MCK_INF, x);
		}
		return ((y % 2) == 0) ? 0.0 : mc_copysign(0.0, x);
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
	return mc_pow(x, mc_cast(double, y));
}

MC_TARGET_FUNC long double mc_powil(const long double x, const int y)
{
	if (x == 0.0L) {
		if (y == 0) {
			return 1.0L;
		} else if (y < 0) {
			return mc_copysignl(MCK_INF, x);
		}
		return ((y % 2) == 0) ? 0.0L : mc_copysignl(0.0L, x);
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
	return mc_powl(x, mc_cast(long double, y));
}

#endif /* !MC_POWI_H */

/* EOF */
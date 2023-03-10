//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rand_uniform_int.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/rand/mc_randu.h>

#ifndef MC_RAND_UNIFORM_INT_H
#define MC_RAND_UNIFORM_INT_H

#pragma mark - mc_rand_uniform_int -

MC_TARGET_FUNC int mc_rand_uniform_int(const int a, const int b)
{
//!# Uniform int distribution generator.
	int r;
	if (mc_iabs(a) < 0x1000001 && mc_iabs(b) < 0x1000001) {
		if (!(a >= b)) {
			const float u = mc_randuf();
			const float c = mc_cast(float, a);
			const float d = mc_cast(float, b);
			r             = mc_cast(int, (c + (u * (d - (c + 1.0f)))));
			r             = r > b ? b : (r < a ? a : r);
		} else {
			r = a;
		}
	} else {
		if (!(a >= b)) {
			const double u = mc_randu();
			const double c = mc_cast(double, a);
			const double d = mc_cast(double, b);
			r              = mc_cast(int, (c + (u * (d - (c + 1.0)))));
			r              = r > b ? b : (r < a ? a : r);
		} else {
			r = a;
		}
	}
	return r;
}

#pragma mark - mc_rand_uniform_int16 -

MC_TARGET_FUNC int16_t mc_rand_uniform_int16(const int16_t a, const int16_t b)
{
//!# Uniform int distribution generator.
	int16_t r;
	if (!(a >= b)) {
		const float u = mc_randuf();
		const float c = mc_cast(float, a);
		const float d = mc_cast(float, b);
		r             = mc_cast(int16_t, (c + (u * (d - (c + 1.0f)))));
		r             = r > b ? b : (r < a ? a : r);
	} else {
		r = a;
	}
	return r;
}

#pragma mark - mc_rand_uniform_int32 -

MC_TARGET_FUNC int32_t mc_rand_uniform_int32(const int32_t a, const int32_t b)
{
//!# Uniform int distribution generator.
	int32_t r;
	if (!(a >= b)) {
		const double u = mc_randu();
		const double c = mc_cast(double, a);
		const double d = mc_cast(double, b);
		r              = mc_cast(int32_t, (c + (u * (d - (c + 1.0)))));
		r              = r > b ? b : (r < a ? a : r);
	} else {
		r = a;
	}
	return r;
}

#pragma mark - mc_rand_uniform_int64 -

MC_TARGET_FUNC int64_t mc_rand_uniform_int64(const int64_t a, const int64_t b)
{
//!# Uniform int distribution generator.
	int64_t r;
	if (!(a >= b)) {
		const long double u = mc_randul();
		const long double c = mc_cast(long double, a);
		const long double d = mc_cast(long double, b);
		r                   = mc_cast(int64_t, (c + (u * (d - (c + 1.0L)))));
		r                   = r > b ? b : (r < a ? a : r);
	} else {
		r = a;
	}
	return r;
}

#endif /* !MC_RAND_UNIFORM_INT_H */

/* EOF */
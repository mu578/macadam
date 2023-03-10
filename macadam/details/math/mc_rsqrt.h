//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_rsqrt.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>

#ifndef MC_RSQRT_H
#define MC_RSQRT_H

#pragma mark - mc_rsqrt -

MC_TARGET_FUNC float mc_rsqrtf(const float x)
{
	if (mc_isnan(x) || mc_isinf(x)) {
		return x;
	}
	if (x <= 0.0f) {
		return MCK_INF;
	}
#	if MC_TARGET_HAVE_SSE
	//__m128 y  = _mm_set_ss(x);
	//__m128 e  = _mm_rsqrt_ss(y);
	//__m128 e3 = _mm_mul_ss(_mm_mul_ss(e, e), e);
	//__m128 h  = _mm_set_ss(0.5f);
	//return _mm_cvtss_f32(_mm_add_ss(e, _mm_mul_ss(h, _mm_sub_ss(e, _mm_mul_ss(y, e3)))));
	const float x2 = x * 0.5f;
	const float f3 = 1.5f;
	float y;
	_mm_store_ss(&y, _mm_rsqrt_ss(_mm_set_ps1(x)));
	return y * (f3 - (x2 * y * y));
#	elif MC_TARGET_HAVE_NEON
	float32x2_t y = vdup_n_f32(x);
	float32x2_t e = vrsqrte_f32(y);
	e             = vmul_f32(e, vrsqrts_f32(y, vmul_f32(e, e)));
	e             = vmul_f32(e, vrsqrts_f32(y, vmul_f32(e, e)));
	return vget_lane_f32(e, 0);
#	elif MC_TARGET_C99 || MC_TARGET_CPP11
	const float x2                      = x * 0.5f;
	const float f3                      = 1.5f;
	union { float f; uint32_t i; } conv = { .f = x };
	//conv.i                              = 0x5F3759DF - (conv.i >> 1);
	conv.i                              = 0x5F375A86 - (conv.i >> 1);
	conv.f                              = conv.f * (f3 - (x2 * conv.f * conv.f));
	return conv.f;
#	else
	return 1.0f / x;
#	endif
} 

MC_TARGET_FUNC double mc_rsqrt(const double x)
{
	if (mc_isnan(x) || mc_isinf(x)) {
		return x;
	}
	if (x <= 0.0) {
		return MCK_INF;
	}
#	if MC_TARGET_HAVE_SSE
		double r;
		const __m128d d = _mm_set_pd(0.0, 1.0);
		__m128d v       = _mm_set_pd(0.0, x);
		v               = _mm_sqrt_pd(v);
		v               = _mm_div_pd(d, v);
		_mm_storel_pd(&r, v);
		return r;
#	elif MC_TARGET_C99 || MC_TARGET_CPP11
	const double x2                      = x * 0.5;
	const double f3                      = 1.5;
	union { double d; uint64_t i; } conv = { .d = x };
	conv.i                               = 0x5FE6EB50C7B537A9 - (conv.i >> 1);
	conv.d                               = conv.d * (f3 - (x2 * conv.d * conv.d));
	return conv.d;
#	else
	return 1.0 / x;
#	endif
}

MC_TARGET_FUNC long double mc_rsqrtl(const long double x)
{
	if (mc_isnan(x) || mc_isinf(x)) {
		return x;
	}
	if (x <= 0.0L) {
		return MCK_INF;
	}
	return 1.0L / x;
}

#endif /* !MC_RSQRT_H */

/* EOF */
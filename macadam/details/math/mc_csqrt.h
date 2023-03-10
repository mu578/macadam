//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_csqrt.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zsqrt.h>

#ifndef MC_CSQRT_H
#define MC_CSQRT_H

#pragma mark - mc_csqrt -

MC_TARGET_PROC mc_complex_float_t mc_csqrtf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return csqrtf(c);
#	else
	mc_complex_float_t z;
	mc_zsqrtf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_csqrt(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return csqrt(c);
#	else
	mc_complex_double_t z;
	mc_zsqrt(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_csqrtl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return csqrtl(c);
#	else
	mc_complex_long_double_t z;
	mc_zsqrtl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CSQRT_H */

/* EOF */
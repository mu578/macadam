//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cdiv.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zdiv.h>

#ifndef MC_CDIV_H
#define MC_CDIV_H

#pragma mark - mc_cdiv -

MC_TARGET_PROC mc_complex_float_t mc_cdivf(const mc_complex_float_t a, const mc_complex_float_t b)
{
#	if MC_TARGET_C99_COMPLEX
		return a / b;
#	else
	mc_complex_float_t c;
	mc_zdivf(&mc_cmplxrf(c), &mc_cmplxif(c)
		, mc_cmplxrf(a), mc_cmplxif(a)
		, mc_cmplxrf(b), mc_cmplxif(b)
	);
	return c;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_cdiv(const mc_complex_double_t a, const mc_complex_double_t b)
{
#	if MC_TARGET_C99_COMPLEX
		return a / b;
#	else
	mc_complex_double_t c;
	mc_zdiv(&mc_cmplxr(c), &mc_cmplxi(c)
		, mc_cmplxr(a), mc_cmplxi(a)
		, mc_cmplxr(b), mc_cmplxi(b)
	);
	return c;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_cdivl(const mc_complex_long_double_t a, const mc_complex_long_double_t b)
{
#	if MC_TARGET_C99_COMPLEX
		return a / b;
#	else
	mc_complex_long_double_t c;
	mc_zdivl(&mc_cmplxrl(c), &mc_cmplxil(c)
		, mc_cmplxrl(a), mc_cmplxil(a)
		, mc_cmplxrl(b), mc_cmplxil(b)
	);
	return c;
#	endif
}

#endif /* !MC_CDIV_H */

/* EOF */
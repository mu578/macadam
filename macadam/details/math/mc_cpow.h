//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cpow.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zpow.h>

#ifndef MC_CPOW_H
#define MC_CPOW_H

#pragma mark - mc_cpow -

MC_TARGET_PROC mc_complex_float_t mc_cpowf(const mc_complex_float_t x, const mc_complex_float_t y)
{
#	if MC_TARGET_C99_COMPLEX
		return cpowf(x, y);
#	else
	mc_complex_float_t c;
	mc_zpowf(&mc_cmplxrf(c), &mc_cmplxif(c)
		, mc_cmplxrf(x), mc_cmplxif(x)
		, mc_cmplxrf(y), mc_cmplxif(y)
	);
	return c;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_cpow(const mc_complex_double_t x, const mc_complex_double_t y)
{
#	if MC_TARGET_C99_COMPLEX
		return cpow(x, y);
#	else
	mc_complex_double_t c;
	mc_zpow(&mc_cmplxr(c), &mc_cmplxi(c)
		, mc_cmplxr(x), mc_cmplxi(x)
		, mc_cmplxr(y), mc_cmplxi(y)
	);
	return c;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_cpowl(const mc_complex_long_double_t x, const mc_complex_long_double_t y)
{
#	if MC_TARGET_C99_COMPLEX
		return cpowl(x, y);
#	else
	mc_complex_long_double_t c;
	mc_zpowl(&mc_cmplxrl(c), &mc_cmplxil(c)
		, mc_cmplxrl(x), mc_cmplxil(x)
		, mc_cmplxrl(y), mc_cmplxil(y)
	);
	return c;
#	endif
}

#endif /* !MC_CPOW_H */

/* EOF */
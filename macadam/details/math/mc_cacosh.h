//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cacosh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zacosh.h>

#ifndef MC_CACOSH_H
#define MC_CACOSH_H

#pragma mark - mc_cacosh -

MC_TARGET_PROC mc_complex_float_t mc_cacoshf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cacoshf(c);
#	else
	mc_complex_float_t z;
	mc_zacoshf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_cacosh(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cacosh(c);
#	else
	mc_complex_double_t z;
	mc_zacosh(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_cacoshl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cacoshl(c);
#	else
	mc_complex_long_double_t z;
	mc_zacoshl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CACOSH_H */

/* EOF */
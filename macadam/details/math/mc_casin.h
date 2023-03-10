//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_casin.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zasin.h>

#ifndef MC_CASIN_H
#define MC_CASIN_H

#pragma mark - mc_casin -

MC_TARGET_PROC mc_complex_float_t mc_casinf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return casinf(c);
#	else
	mc_complex_float_t z;
	mc_zasinf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_casin(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return casin(c);
#	else
	mc_complex_double_t z;
	mc_zasin(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_casinl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return casinl(c);
#	else
	mc_complex_long_double_t z;
	mc_zasinl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CASIN_H */

/* EOF */
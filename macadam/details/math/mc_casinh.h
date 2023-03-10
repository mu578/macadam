//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_casinh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zasinh.h>

#ifndef MC_CASINH_H
#define MC_CASINH_H

#pragma mark - mc_casinh -

MC_TARGET_PROC mc_complex_float_t mc_casinhf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return casinhf(c);
#	else
	mc_complex_float_t z;
	mc_zasinhf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_casinh(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return casinh(c);
#	else
	mc_complex_double_t z;
	mc_zasinh(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_casinhl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return casinhl(c);
#	else
	mc_complex_long_double_t z;
	mc_zasinhl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CASINH_H */

/* EOF */
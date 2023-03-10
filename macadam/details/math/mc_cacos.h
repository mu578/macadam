//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cacos.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zacos.h>

#ifndef MC_CACOS_H
#define MC_CACOS_H

#pragma mark - mc_cacos -

MC_TARGET_PROC mc_complex_float_t mc_cacosf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cacosf(c);
#	else
	mc_complex_float_t z;
	mc_zacosf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_cacos(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cacos(c);
#	else
	mc_complex_double_t z;
	mc_zacos(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_cacosl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cacosl(c);
#	else
	mc_complex_long_double_t z;
	mc_zacosl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CACOS_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_catanh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zatanh.h>

#ifndef MC_CATANH_H
#define MC_CATANH_H

#pragma mark - mc_catanh -

MC_TARGET_PROC mc_complex_float_t mc_catanhf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return catanhf(c);
#	else
	mc_complex_float_t z;
	mc_zatanhf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_catanh(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return catanh(c);
#	else
	mc_complex_double_t z;
	mc_zatanh(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_catanhl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return catanhl(c);
#	else
	mc_complex_long_double_t z;
	mc_zatanhl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CATANH_H */

/* EOF */
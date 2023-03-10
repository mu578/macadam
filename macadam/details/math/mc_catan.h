//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_catan.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zatan.h>

#ifndef MC_CATAN_H
#define MC_CATAN_H

#pragma mark - mc_catan -

MC_TARGET_PROC mc_complex_float_t mc_catanf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return catanf(c);
#	else
	mc_complex_float_t z;
	mc_zatanf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_catan(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return catan(c);
#	else
	mc_complex_double_t z;
	mc_zatan(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_catanl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return catanl(c);
#	else
	mc_complex_long_double_t z;
	mc_zatanl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CATAN_H */

/* EOF */
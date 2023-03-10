//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_csin.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zsin.h>

#ifndef MC_CSIN_H
#define MC_CSIN_H

#pragma mark - mc_csin -

MC_TARGET_PROC mc_complex_float_t mc_csinf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return csinf(c);
#	else
	mc_complex_float_t z;
	mc_zsinf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_csin(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return csin(c);
#	else
	mc_complex_double_t z;
	mc_zsin(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_csinl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return csinl(c);
#	else
	mc_complex_long_double_t z;
	mc_zsinl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CSIN_H */

/* EOF */
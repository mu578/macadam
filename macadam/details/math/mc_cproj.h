//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cproj.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zproj.h>

#ifndef MC_CPROJ_H
#define MC_CPROJ_H

#pragma mark - mc_cproj -

MC_TARGET_PROC mc_complex_float_t mc_cprojf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cprojf(c);
#	else
	mc_complex_float_t z;
	mc_zprojf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_cproj(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cproj(c);
#	else
	mc_complex_double_t z;
	mc_zproj(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_cprojl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cprojl(c);
#	else
	mc_complex_long_double_t z;
	mc_zprojl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CPROJ_H */

/* EOF */
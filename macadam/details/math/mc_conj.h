//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_conj.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_CONJ_H
#define MC_CONJ_H

#pragma mark - mc_conj -

MC_TARGET_PROC mc_complex_float_t mc_conjf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return conjf(c);
#	else
	mc_complex_float_t z;
	mc_zconjf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_conj(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return conj(c);
#	else
	mc_complex_double_t z;
	mc_zconj(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_conjl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return conjl(c);
#	else
	mc_complex_long_double_t z;
	mc_zconjl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CONJ_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cabs.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zabs.h>

#ifndef MC_CABS_H
#define MC_CABS_H

#pragma mark - mc_cabs -

MC_TARGET_PROC float mc_cabsf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cabsf(c);
#	else
		return mc_zabsf(mc_cmplxrf(c), mc_cmplxif(c));
#	endif
}

MC_TARGET_PROC double mc_cabs(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cabs(c);
#	else
		return mc_zabs(mc_cmplxr(c), mc_cmplxi(c));
#	endif
}

MC_TARGET_PROC long double mc_cabsl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cabsl(c);
#	else
		return mc_zabsl(mc_cmplxrl(c), mc_cmplxil(c));
#	endif
}

#endif /* !MC_CABS_H */

/* EOF */
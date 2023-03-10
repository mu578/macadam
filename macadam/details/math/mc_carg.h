//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_carg.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zarg.h>

#ifndef MC_CARG_H
#define MC_CARG_H

#pragma mark - mc_carg -

MC_TARGET_PROC float mc_cargf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cargf(c);
#	else
	return mc_zargf(mc_cmplxrf(c), mc_cmplxif(c));
#	endif
}

MC_TARGET_PROC double mc_carg(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return carg(c);
#	else
	return mc_zarg(mc_cmplxr(c), mc_cmplxi(c));
#	endif
}

MC_TARGET_PROC long double mc_cargl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cargl(c);
#	else
	return mc_zargl(mc_cmplxrl(c), mc_cmplxil(c));
#	endif
}

#endif /* !MC_CARG_H */

/* EOF */
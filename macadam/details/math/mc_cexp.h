//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cexp.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zexp.h>

#ifndef MC_CEXP_H
#define MC_CEXP_H

#pragma mark - mc_cexp -

MC_TARGET_PROC mc_complex_float_t mc_cexpf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cexpf(c);
#	else
	mc_complex_float_t z;
	mc_zexpf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_cexp(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cexp(c);
#	else
	mc_complex_double_t z;
	mc_zexp(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_cexpl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return cexpl(c);
#	else
	mc_complex_long_double_t z;
	mc_zexpl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CEXP_H */

/* EOF */
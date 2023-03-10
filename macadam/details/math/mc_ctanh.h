//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ctanh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_ztanh.h>

#ifndef MC_CTANH_H
#define MC_CTANH_H

#pragma mark - mc_ctanh -

MC_TARGET_PROC mc_complex_float_t mc_ctanhf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return ctanhf(c);
#	else
	mc_complex_float_t z;
	mc_ztanhf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_ctanh(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return ctanh(c);
#	else
	mc_complex_double_t z;
	mc_ztanh(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_ctanhl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return ctanhl(c);
#	else
	mc_complex_long_double_t z;
	mc_ztanhl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CTANH_H */

/* EOF */
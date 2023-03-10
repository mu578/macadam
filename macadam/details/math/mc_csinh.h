//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_csinh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zsinh.h>

#ifndef MC_CSINH_H
#define MC_CSINH_H

#pragma mark - mc_csinh -

MC_TARGET_PROC mc_complex_float_t mc_csinhf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return csinhf(c);
#	else
	mc_complex_float_t z;
	mc_zsinhf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_csinh(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return csinh(c);
#	else
	mc_complex_double_t z;
	mc_zsinh(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_csinhl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return csinhl(c);
#	else
	mc_complex_long_double_t z;
	mc_zsinhl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CSINH_H */

/* EOF */
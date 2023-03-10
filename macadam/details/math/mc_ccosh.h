//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ccosh.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zcosh.h>

#ifndef MC_CCOSH_H
#define MC_CCOSH_H

#pragma mark - mc_ccosh -

MC_TARGET_PROC mc_complex_float_t mc_ccoshf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return ccoshf(c);
#	else
	mc_complex_float_t z;
	mc_zcoshf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_ccosh(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return ccosh(c);
#	else
	mc_complex_double_t z;
	mc_zcosh(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_ccoshl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return ccoshl(c);
#	else
	mc_complex_long_double_t z;
	mc_zcoshl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CCOSH_H */

/* EOF */
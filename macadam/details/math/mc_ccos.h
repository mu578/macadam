//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ccos.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zcos.h>

#ifndef MC_CCOS_H
#define MC_CCOS_H

#pragma mark - mc_ccos -

MC_TARGET_PROC mc_complex_float_t mc_ccosf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return ccosf(c);
#	else
	mc_complex_float_t z;
	mc_zcosf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_ccos(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return ccos(c);
#	else
	mc_complex_double_t z;
	mc_zcos(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_ccosl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return ccosl(c);
#	else
	mc_complex_long_double_t z;
	mc_zcosl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CCOS_H */

/* EOF */
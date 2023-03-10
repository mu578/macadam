//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ctan.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_ztan.h>

#ifndef MC_CTAN_H
#define MC_CTAN_H

#pragma mark - mc_ctan -

MC_TARGET_PROC mc_complex_float_t mc_ctanf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return ctanf(c);
#	else
	mc_complex_float_t z;
	mc_ztanf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_ctan(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return ctan(c);
#	else
	mc_complex_double_t z;
	mc_ztan(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_ctanl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return ctanl(c);
#	else
	mc_complex_long_double_t z;
	mc_ztanl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CTAN_H */

/* EOF */
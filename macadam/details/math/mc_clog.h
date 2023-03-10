//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_clog.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zlog.h>

#ifndef MC_CLOG_H
#define MC_CLOG_H

#pragma mark - mc_clog -

MC_TARGET_PROC mc_complex_float_t mc_clogf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return clogf(c);
#	else
	mc_complex_float_t z;
	mc_zlogf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_clog(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return clog(c);
#	else
	mc_complex_double_t z;
	mc_zlog(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_clogl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
		return clogl(c);
#	else
	mc_complex_long_double_t z;
	mc_zlogl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CLOG_H */

/* EOF */
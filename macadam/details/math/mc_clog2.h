//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_clog2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zlog2.h>

#ifndef MC_CLOG2_H
#define MC_CLOG2_H

#pragma mark - mc_clog2 -

MC_TARGET_PROC mc_complex_float_t mc_clog2f(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX && ((defined(__linux__) || defined(__gnu_linux__)) && defined(_GNU_SOURCE))
		return clog2f(c);
#	elif MC_TARGET_C99_COMPLEX
	float z_r, z_i;
	mc_zlog2f(&z_r, &z_i
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return mc_cmplxf(z_r, z_i);
#	else
	mc_complex_float_t z;
	mc_zlog2f(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_clog2(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX && ((defined(__linux__) || defined(__gnu_linux__)) && defined(_GNU_SOURCE))
		return clog2(c);
#	elif MC_TARGET_C99_COMPLEX
	double z_r, z_i;
	mc_zlog2(&z_r, &z_i
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return mc_cmplx(z_r, z_i);
#	else
	mc_complex_double_t z;
	mc_zlog2(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_clog2l(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX && ((defined(__linux__) || defined(__gnu_linux__)) && defined(_GNU_SOURCE))
		return clog2l(c);
#	elif MC_TARGET_C99_COMPLEX
	long double z_r, z_i;
	mc_zlog2l(&z_r, &z_i
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return mc_cmplxl(z_r, z_i);
#	else
	mc_complex_long_double_t z;
	mc_zlog2l(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CLOG2_H */

/* EOF */
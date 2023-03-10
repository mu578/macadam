//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_clog1p.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zlog1p.h>

#ifndef MC_CLOG1P_H
#define MC_CLOG1P_H

#pragma mark - mc_clog1p -

MC_TARGET_PROC mc_complex_float_t mc_clog1pf(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
	float z_r, z_i;
	mc_zlog1pf(&z_r, &z_i
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return mc_cmplxf(z_r, z_i);
#	else
	mc_complex_float_t z;
	mc_zlog1pf(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_clog1p(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
	double z_r, z_i;
	mc_zlog1p(&z_r, &z_i
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return mc_cmplx(z_r, z_i);
#	else
	mc_complex_double_t z;
	mc_zlog1p(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_clog1pl(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
	long double z_r, z_i;
	mc_zlog1pl(&z_r, &z_i
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return mc_cmplxl(z_r, z_i);
#	else
	mc_complex_long_double_t z;
	mc_zlog1pl(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CLOG1P_H */

/* EOF */
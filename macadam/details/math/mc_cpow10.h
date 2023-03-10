//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cpow10.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zpow10.h>

#ifndef MC_CPOW10_H
#define MC_CPOW10_H

#pragma mark - mc_cpow10 -

MC_TARGET_PROC mc_complex_float_t mc_cpow10f(const mc_complex_float_t c)
{
#	if MC_TARGET_C99_COMPLEX
	float z_r, z_i;
	mc_zpow10f(&z_r, &z_i
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return mc_cmplxf(z_r, z_i);
#	else
	mc_complex_float_t z;
	mc_zpow10f(&mc_cmplxrf(z), &mc_cmplxif(z)
		, mc_cmplxrf(c), mc_cmplxif(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_double_t mc_cpow10(const mc_complex_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
	double z_r, z_i;
	mc_zpow10(&z_r, &z_i
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return mc_cmplx(z_r, z_i);
#	else
	mc_complex_double_t z;
	mc_zpow10(&mc_cmplxr(z), &mc_cmplxi(z)
		, mc_cmplxr(c), mc_cmplxi(c)
	);
	return z;
#	endif
}

MC_TARGET_PROC mc_complex_long_double_t mc_cpow10l(const mc_complex_long_double_t c)
{
#	if MC_TARGET_C99_COMPLEX
	long double z_r, z_i;
	mc_zpow10l(&z_r, &z_i
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return mc_cmplxl(z_r, z_i);
#	else
	mc_complex_long_double_t z;
	mc_zpow10l(&mc_cmplxrl(z), &mc_cmplxil(z)
		, mc_cmplxrl(c), mc_cmplxil(c)
	);
	return z;
#	endif
}

#endif /* !MC_CPOW10_H */

/* EOF */
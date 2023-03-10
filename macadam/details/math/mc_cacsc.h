//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cacsc.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zacsc.h>

#ifndef MC_CACSC_H
#define MC_CACSC_H

#pragma mark - mc_cacsc -

MC_TARGET_FUNC mc_complex_float_t mc_cacscf(const mc_complex_float_t c)
{
	mc_complex_float_t z;
	float zr, zi;

	mc_zacscf(&zr, &zi, mc_cmplxrf(c), mc_cmplxif(c));
	z = mc_cmplxf(zr, zi);

	return z;
}

MC_TARGET_FUNC mc_complex_double_t mc_cacsc(const mc_complex_double_t c)
{
	mc_complex_double_t z;
	double zr, zi;

	mc_zacsc(&zr, &zi, mc_cmplxr(c), mc_cmplxi(c));
	z = mc_cmplx(zr, zi);

	return z;
}

MC_TARGET_FUNC mc_complex_long_double_t mc_cacscl(const mc_complex_long_double_t c)
{
	mc_complex_long_double_t z;
	long double zr, zi;

	mc_zacscl(&zr, &zi, mc_cmplxrl(c), mc_cmplxil(c));
	z = mc_cmplxl(zr, zi);

	return z;
}

#endif /* !MC_CACSC_H */

/* EOF */
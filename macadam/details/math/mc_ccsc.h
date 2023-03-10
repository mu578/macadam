//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ccsc.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zcsc.h>

#ifndef MC_CCSC_H
#define MC_CCSC_H

#pragma mark - mc_ccsc -

MC_TARGET_FUNC mc_complex_float_t mc_ccscf(const mc_complex_float_t c)
{
	mc_complex_float_t z;
	float zr, zi;

	mc_zcscf(&zr, &zi, mc_cmplxrf(c), mc_cmplxif(c));
	z = mc_cmplxf(zr, zi);

	return z;
}

MC_TARGET_FUNC mc_complex_double_t mc_ccsc(const mc_complex_double_t c)
{
	mc_complex_double_t z;
	double zr, zi;

	mc_zcsc(&zr, &zi, mc_cmplxr(c), mc_cmplxi(c));
	z = mc_cmplx(zr, zi);

	return z;
}

MC_TARGET_FUNC mc_complex_long_double_t mc_ccscl(const mc_complex_long_double_t c)
{
	mc_complex_long_double_t z;
	long double zr, zi;

	mc_zcscl(&zr, &zi, mc_cmplxrl(c), mc_cmplxil(c));
	z = mc_cmplxl(zr, zi);

	return z;
}

#endif /* !MC_CCSC_H */

/* EOF */
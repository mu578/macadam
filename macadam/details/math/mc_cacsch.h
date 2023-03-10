//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cacsch.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zacsch.h>

#ifndef MC_CACSCH_H
#define MC_CACSCH_H

#pragma mark - mc_cacsch -

MC_TARGET_FUNC mc_complex_float_t mc_cacschf(const mc_complex_float_t c)
{
	mc_complex_float_t z;
	float zr, zi;

	mc_zacschf(&zr, &zi, mc_cmplxrf(c), mc_cmplxif(c));
	z = mc_cmplxf(zr, zi);

	return z;
}

MC_TARGET_FUNC mc_complex_double_t mc_cacsch(const mc_complex_double_t c)
{
	mc_complex_double_t z;
	double zr, zi;

	mc_zacsch(&zr, &zi, mc_cmplxr(c), mc_cmplxi(c));
	z = mc_cmplx(zr, zi);

	return z;
}

MC_TARGET_FUNC mc_complex_long_double_t mc_cacschl(const mc_complex_long_double_t c)
{
	mc_complex_long_double_t z;
	long double zr, zi;

	mc_zacschl(&zr, &zi, mc_cmplxrl(c), mc_cmplxil(c));
	z = mc_cmplxl(zr, zi);

	return z;
}

#endif /* !MC_CACSCH_H */

/* EOF */
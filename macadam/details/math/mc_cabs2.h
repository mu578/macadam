//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cabs2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zabs2.h>

#ifndef MC_CABS2_H
#define MC_CABS2_H

#pragma mark - mc_cabs2 -

MC_TARGET_PROC float mc_cabs2f(const mc_complex_float_t c)
{
	return mc_zabs2f(mc_cmplxrf(c), mc_cmplxif(c));
}

MC_TARGET_PROC double mc_cabs2(const mc_complex_double_t c)
{
	return mc_zabs2(mc_cmplxr(c), mc_cmplxi(c));
}

MC_TARGET_PROC long double mc_cabs2l(const mc_complex_long_double_t c)
{
	return mc_zabs2l(mc_cmplxrl(c), mc_cmplxil(c));
}

#endif /* !MC_CABS2_H */

/* EOF */
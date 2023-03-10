//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cimag.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zmul.h>

#ifndef MC_CIMAG_H
#define MC_CIMAG_H

#pragma mark - mc_cimag -

MC_TARGET_PROC float mc_cimagf(const mc_complex_float_t c)
{
	return mc_cmplxif(c);
}

MC_TARGET_PROC double mc_cimag(const mc_complex_double_t c)
{
	return mc_cmplxi(c);
}

MC_TARGET_PROC long double mc_cimagl(const mc_complex_long_double_t c)
{
	return mc_cmplxil(c);
}

#endif /* !MC_CIMAG_H */

/* EOF */
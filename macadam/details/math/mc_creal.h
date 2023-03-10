//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_creal.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_zmul.h>

#ifndef MC_CREAL_H
#define MC_CREAL_H

#pragma mark - mc_creal -

MC_TARGET_PROC float mc_crealf(const mc_complex_float_t c)
{
	return mc_cmplxrf(c);
}

MC_TARGET_PROC double mc_creal(const mc_complex_double_t c)
{
	return mc_cmplxr(c);
}

MC_TARGET_PROC long double mc_creall(const mc_complex_long_double_t c)
{
	return mc_cmplxrl(c);
}

#endif /* !MC_CREAL_H */

/* EOF */
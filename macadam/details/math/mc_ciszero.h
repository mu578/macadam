//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ciszero.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_CISZERO_H
#define MC_CISZERO_H

#pragma mark - mc_ciszero -

MC_TARGET_PROC int mc_ciszerof(const mc_complex_float_t c)
{
	return (mc_cmplxrf(c) == 0.0f && mc_cmplxif(c) == 0.0f) ? 1 : 0;
}

MC_TARGET_PROC int mc_ciszero(const mc_complex_double_t c)
{
	return (mc_cmplxr(c) == 0.0 && mc_cmplxi(c) == 0.0) ? 1 : 0;
}

MC_TARGET_PROC int mc_ciszerol(const mc_complex_long_double_t c)
{
	return (mc_cmplxrl(c) == 0.0L && mc_cmplxil(c) == 0.0L) ? 1 : 0;
}

#endif /* !MC_CISZERO_H */

/* EOF */
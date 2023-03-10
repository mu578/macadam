//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ciseq.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_CISEQ_H
#define MC_CISEQ_H

#pragma mark - mc_ciseq -

MC_TARGET_PROC int mc_ciseqf(const mc_complex_float_t a, const mc_complex_float_t b)
{
	return (mc_cmplxrf(a) == mc_cmplxrf(b) && mc_cmplxif(a) == mc_cmplxif(b)) ? 1 : 0;
}

MC_TARGET_PROC int mc_ciseq(const mc_complex_double_t a, const mc_complex_double_t b)
{
	return (mc_cmplxr(a) == mc_cmplxr(b) && mc_cmplxi(a) == mc_cmplxi(b)) ? 1 : 0;
}

MC_TARGET_PROC int mc_ciseql(const mc_complex_long_double_t a, const mc_complex_long_double_t b)
{
	return (mc_cmplxrl(a) == mc_cmplxrl(b) && mc_cmplxil(a) == mc_cmplxil(b)) ? 1 : 0;
}

#endif /* !MC_CISEQ_H */

/* EOF */
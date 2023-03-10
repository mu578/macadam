//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ziseq.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ZISEQ_H
#define MC_ZISEQ_H

#pragma mark - mc_ziseq -

MC_TARGET_PROC int mc_ziseqf(const float a_r, const float a_i, const float b_r, const float b_i)
{
	return (a_r == b_r && a_i == b_i) ? 1 : 0;
}

MC_TARGET_PROC int mc_ziseq(const double a_r, double a_i, const double b_r, const double b_i)
{
	return (a_r == b_r && a_i == b_i) ? 1 : 0;
}

MC_TARGET_PROC int mc_ziseql(const long double a_r, const long double a_i, const long double b_r, const long double b_i)
{
	return (a_r == b_r && a_i == b_i) ? 1 : 0;
}

#endif /* !MC_ZISEQ_H */

/* EOF */
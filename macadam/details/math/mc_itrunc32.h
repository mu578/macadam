//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_itrunc32.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ITRUNC32_H
#define MC_ITRUNC32_H

#pragma mark - mc_itrunc32 -

MC_TARGET_PROC int32_t mc_itrunc32f(const float x)
{
	return mc_cast(int32_t, x);
}

MC_TARGET_PROC int32_t mc_itrunc32(const double x)
{
	return mc_cast(int32_t, x);
}

MC_TARGET_PROC int32_t mc_itrunc32l(const long double x)
{
	return mc_cast(int32_t, x);
}

#endif /* !MC_ITRUNC32_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_itrunc16.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ITRUNC16_H
#define MC_ITRUNC16_H

#pragma mark - mc_itrunc16 -

MC_TARGET_PROC int16_t mc_itrunc16f(const float x)
{
	return mc_cast(int16_t, x);
}

MC_TARGET_PROC int16_t mc_itrunc16(const double x)
{
	return mc_cast(int16_t, x);
}

MC_TARGET_PROC int16_t mc_itrunc16l(const long double x)
{
	return mc_cast(int16_t, x);
}

#endif /* !MC_ITRUNC16_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ziszero.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ZISZERO_H
#define MC_ZISZERO_H

#pragma mark - mc_ziszero -

MC_TARGET_PROC int mc_ziszerof(const float z_r, const float z_i)
{
	return (z_r == 0.0f && z_i == 0.0f) ? 1 : 0;
}

MC_TARGET_PROC int mc_ziszero(const double z_r, const double z_i)
{
	return (z_r == 0.0 && z_i == 0.0) ? 1 : 0;
}

MC_TARGET_PROC int mc_ziszerol(const long double z_r, const long double z_i)
{
	return (z_r == 0.0L && z_i == 0.0L) ? 1 : 0;
}

#endif /* !MC_ZISZERO_H */

/* EOF */
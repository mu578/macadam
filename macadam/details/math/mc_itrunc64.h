//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_itrunc64.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ITRUNC64_H
#define MC_ITRUNC64_H

#pragma mark - mc_itrunc64 -

MC_TARGET_PROC int64_t mc_itrunc64f(const float x)
{
	return mc_cast(int64_t, x);
}

MC_TARGET_PROC int64_t mc_itrunc64(const double x)
{
	return mc_cast(int64_t, x);
}

MC_TARGET_PROC int64_t mc_itrunc64l(const long double x)
{
	return mc_cast(int64_t, x);
}

#endif /* !MC_ITRUNC64_H */

/* EOF */
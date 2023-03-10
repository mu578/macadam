//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_raise2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_RAISE2_H
#define MC_RAISE2_H

#pragma mark - mc_raise2 -

MC_TARGET_PROC float mc_raise2f(const float x)
{
	return x * x;
}

MC_TARGET_PROC double mc_raise2(const double x)
{
	return x * x;
}

MC_TARGET_PROC long double mc_raise2l(const long double x)
{
	return x * x;
}

#endif /* !MC_RAISE2_H */

/* EOF */
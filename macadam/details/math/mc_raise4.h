//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_raise4.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_RAISE4_H
#define MC_RAISE4_H

#pragma mark - mc_raise4 -

MC_TARGET_PROC float mc_raise4f(const float x)
{
	return x * x * x * x;
}

MC_TARGET_PROC double mc_raise4(const double x)
{
	return x * x * x * x;
}

MC_TARGET_PROC long double mc_raise4l(const long double x)
{
	return x * x * x * x;
}

#endif /* !MC_RAISE4_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_raise5.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_RAISE5_H
#define MC_RAISE5_H

#pragma mark - mc_raise5 -

MC_TARGET_PROC float mc_raise5f(const float x)
{
	return x * x * x * x * x;
}

MC_TARGET_PROC double mc_raise5(const double x)
{
	return x * x * x * x * x;
}

MC_TARGET_PROC long double mc_raise5l(const long double x)
{
	return x * x * x * x * x;
}

#endif /* !MC_RAISE5_H */

/* EOF */
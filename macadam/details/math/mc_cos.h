//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_cos.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_COS_H
#define MC_COS_H

#pragma mark - mc_cos -

MC_TARGET_FUNC float mc_cosf(const float x)
{
#	if MC_TARGET_CPP98
	return ::cosf(x);
#	else
	return cosf(x);
#	endif
}

MC_TARGET_FUNC double mc_cos(const double x)
{
#	if MC_TARGET_CPP98
	return ::cos(x);
#	else
	return cos(x);
#	endif
}

MC_TARGET_FUNC long double mc_cosl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::cosl(x);
#	else
	return cosl(x);
#	endif
}

#endif /* !MC_COS_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_scalbn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef SCALBN_H
#define SCALBN_H

#pragma mark - mc_scalbn -

MC_TARGET_FUNC float mc_scalbnf(const float x, const int y)
{
#	if MC_TARGET_CPP98
	return ::scalbnf(x, y);
#	else
	return scalbnf(x, y);
#	endif
}

MC_TARGET_FUNC double mc_scalbn(const double x, const int y)
{
#	if MC_TARGET_CPP98
	return ::scalbn(x, y);
#	else
	return scalbn(x, y);
#	endif
}

MC_TARGET_FUNC long double mc_scalbnl(const long double x, const int y)
{
#	if MC_TARGET_CPP98
	return ::scalbnl(x, y);
#	else
	return scalbnl(x, y);
#	endif
}

#endif /* !SCALBN_H */

/* EOF */
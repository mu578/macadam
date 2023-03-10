//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_modf.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_MODF_H
#define MC_MODF_H

#pragma mark - mc_modf -

MC_TARGET_FUNC float mc_modff(const float x, float * y)
{
#	if MC_TARGET_CPP98
	return ::modff(x, y);
#	else
	return modff(x, y);
#	endif
}

MC_TARGET_FUNC double mc_modf(const double x, double * y)
{
#	if MC_TARGET_CPP98
	return ::modf(x, y);
#	else
	return modf(x, y);
#	endif
}

MC_TARGET_FUNC long double mc_modfl(const long double x, long double * y)
{
#	if MC_TARGET_CPP98
	return ::modfl(x, y);
#	else
	return modfl(x, y);
#	endif
}

#endif /* !MC_MODF_H */

/* EOF */
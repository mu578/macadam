//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_erfc.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ERFC_H
#define MC_ERFC_H

#pragma mark - mc_erfc -

MC_TARGET_FUNC float mc_erfcf(const float x)
{
#	if MC_TARGET_CPP98
	return ::erfcf(x);
#	else
	return erfcf(x);
#	endif
}

MC_TARGET_FUNC double mc_erfc(const double x)
{
#	if MC_TARGET_CPP98
	return ::erfc(x);
#	else
	return erfc(x);
#	endif
}

MC_TARGET_FUNC long double mc_erfcl(const long double x)
{
#	if MC_TARGET_CPP98
	return ::erfcl(x);
#	else
	return erfcl(x);
#	endif
}

#endif /* !MC_ERFC_H */

/* EOF */
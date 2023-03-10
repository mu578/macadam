//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_isnan.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ISNAN_H
#define MC_ISNAN_H

#pragma mark - mc_isnan -

#	if MC_TARGET_CPP98
#	define mc_isnan(x) (::isnan(x) != 0 ? 1 : 0)
#	else
#	define mc_isnan(x) (isnan(x) != 0 ? 1 : 0)
#	endif

#endif /* !MC_ISNAN_H */

/* EOF */
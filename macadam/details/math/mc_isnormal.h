//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_isnormal.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ISNORMAL_H
#define MC_ISNORMAL_H

#pragma mark - mc_isnormal -

#	if MC_TARGET_CPP98
#	define mc_isnormal(x) (::isnormal(x) != 0 ? 1 : 0)
#	else
#	define mc_isnormal(x) (isnormal(x) != 0 ? 1 : 0)
#	endif

#endif /* !MC_ISNORMAL_H */

/* EOF */
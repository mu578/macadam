//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_isfinite.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ISFINITE_H
#define MC_ISFINITE_H

#pragma mark - mc_isfinite -

#	if MC_TARGET_CPP98
#	define mc_isfinite(x) (::isfinite(x) ? 1 : 0)
#	else
#	define mc_isfinite(x) (isfinite(x) ? 1 : 0)
#	endif

#endif /* !MC_ISFINITE_H */

/* EOF */
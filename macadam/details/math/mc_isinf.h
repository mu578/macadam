//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_isinf.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ISINF_H
#define MC_ISINF_H

#pragma mark - mc_isinf -

#	if MC_TARGET_CPP98
#	define mc_isinf(x) (::isinf(x) != 0 ? 1 : 0)
#	else
#	define mc_isinf(x) (isinf(x) != 0 ? 1 : 0)
#	endif

#	if MC_TARGET_CPP98
#	define mc_isinfp(x) (::isinf(x) > 0 || x == MCK_INFP ? 1 : 0)
#	else
#	define mc_isinfp(x) (isinf(x) > 0 || x == MCK_INFP ? 1 : 0)
#	endif

#	if MC_TARGET_CPP98
#	define mc_isinfn(x) (::isinf(x) < 0 || x == MCK_INFN ? 1 : 0)
#	else
#	define mc_isinfn(x) (isinf(x) < 0 || x == MCK_INFN ? 1 : 0)
#	endif

#endif /* !MC_ISINF_H */

/* EOF */
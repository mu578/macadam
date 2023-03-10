//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_signbit.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_SIGNBIT_H
#define MC_SIGNBIT_H

#pragma mark - mc_signbit -

#	define mc_signbitf(x) MC_TARGET_SIGNBITF(x)
#	define mc_signbit(x)  MC_TARGET_SIGNBIT(x)
#	define mc_signbitl(x) MC_TARGET_SIGNBITL(x)

#endif /* !MC_SIGNBIT_H */

/* EOF */
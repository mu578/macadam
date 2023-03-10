//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ones1x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_ONES1X2_H
#define MC_ONES1X2_H

#pragma mark - mc_ones1x2 -

MC_TARGET_FUNC void mc_ones1x2f(float x[2])
{
	x[0] = 1.0f; x[1] = 1.0f;
}

MC_TARGET_FUNC void mc_ones1x2(double x[2])
{
	x[0] = 1.0; x[1] = 1.0;
}

MC_TARGET_FUNC void mc_ones1x2l(long double x[2])
{
	x[0] = 1.0L; x[1] = 1.0L;
}

#endif /* !MC_ONES1X2_H */

/* EOF */
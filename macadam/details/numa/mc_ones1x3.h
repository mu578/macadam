//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ones1x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_ONES1X3_H
#define MC_ONES1X3_H

#pragma mark - mc_ones1x3 -

MC_TARGET_FUNC void mc_ones1x3f(float x[3])
{
	x[0] = 1.0f; x[1] = 1.0f; x[2] = 1.0f;
}

MC_TARGET_FUNC void mc_ones1x3(double x[3])
{
	x[0] = 1.0; x[1] = 1.0; x[2] = 1.0;
}

MC_TARGET_FUNC void mc_ones1x3l(long double x[3])
{
	x[0] = 1.0L; x[1] = 1.0L; x[2] = 1.0L;
}

#endif /* !MC_ONES1X3_H */

/* EOF */
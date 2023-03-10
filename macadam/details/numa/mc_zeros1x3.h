//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zeros1x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_ZEROS1X3_H
#define MC_ZEROS1X3_H

#pragma mark - mc_zeros1x3 -

MC_TARGET_FUNC void mc_zeros1x3f(float x[3])
{
	x[0] = 0.0f; x[1] = 0.0f; x[2] = 0.0f;
}

MC_TARGET_FUNC void mc_zeros1x3(double x[3])
{
	x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
}

MC_TARGET_FUNC void mc_zeros1x3l(long double x[3])
{
	x[0] = 0.0L; x[1] = 0.0L; x[2] = 0.0L;
}

#endif /* !MC_ZEROS1X3_H */

/* EOF */
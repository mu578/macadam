//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_ones2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_ONES2X2_H
#define MC_ONES2X2_H

#pragma mark - mc_ones2x2 -

MC_TARGET_FUNC void mc_ones2x2f(float a[4])
{
	a[0] = 1.0f; a[1] = 1.0f;
	a[2] = 1.0f; a[3] = 1.0f;
}

MC_TARGET_FUNC void mc_ones2x2(double a[4])
{
	a[0] = 1.0; a[1] = 1.0;
	a[2] = 1.0; a[3] = 1.0;
}

MC_TARGET_FUNC void mc_ones2x2l(long double a[4])
{
	a[0] = 1.0L; a[1] = 1.0L;
	a[2] = 1.0L; a[3] = 1.0L;
}

#endif /* !MC_ONES2X2_H */

/* EOF */
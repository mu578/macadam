//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_zeros2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_ZEROS2X2_H
#define MC_ZEROS2X2_H

#pragma mark - mc_zeros2x2 -

MC_TARGET_FUNC void mc_zeros2x2f(float a[4])
{
	a[0] = 0.0f; a[1] = 0.0f;
	a[2] = 0.0f; a[3] = 0.0f;
}

MC_TARGET_FUNC void mc_zeros2x2(double a[4])
{
	a[0] = 0.0; a[1] = 0.0;
	a[2] = 0.0; a[3] = 0.0;
}

MC_TARGET_FUNC void mc_zeros2x2l(long double a[4])
{
	a[0] = 0.0L; a[1] = 0.0L;
	a[2] = 0.0L; a[3] = 0.0L;
}

#endif /* !MC_ZEROS2X2_H */

/* EOF */
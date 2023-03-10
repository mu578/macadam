//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_magic2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_MAGIC2X2_H
#define MC_MAGIC2X2_H

#pragma mark - mc_magic2x2 -

MC_TARGET_FUNC void mc_magic2x2f(float a[4])
{
	a[0] = 4.0f; a[1] = 3.0f;
	a[2] = 1.0f; a[3] = 2.0f;
}

MC_TARGET_FUNC void mc_magic2x2(double a[4])
{
	a[0] = 4.0; a[1] = 3.0;
	a[2] = 1.0; a[3] = 2.0;
}

MC_TARGET_FUNC void mc_magic2x2l(long double a[4])
{
	a[0] = 4.0L; a[1] = 3.0L;
	a[2] = 1.0L; a[3] = 2.0L;
}

#endif /* !MC_MAGIC2X2_H */

/* EOF */
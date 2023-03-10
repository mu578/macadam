//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_magic3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_MAGIC3X3_H
#define MC_MAGIC3X3_H

#pragma mark - mc_magic3x3 -

MC_TARGET_FUNC void mc_magic3x3f(float a[9])
{
	a[0] = 8.0f; a[1] = 1.0f; a[2] = 6.0f;
	a[3] = 3.0f; a[4] = 5.0f; a[5] = 7.0f;
	a[6] = 4.0f; a[7] = 9.0f; a[8] = 2.0f;
}

MC_TARGET_FUNC void mc_magic3x3(double a[9])
{
	a[0] = 8.0; a[1] = 1.0; a[2] = 6.0;
	a[3] = 3.0; a[4] = 5.0; a[5] = 7.0;
	a[6] = 4.0; a[7] = 9.0; a[8] = 2.0;
}

MC_TARGET_FUNC void mc_magic3x3l(long double a[9])
{
	a[0] = 8.0L; a[1] = 1.0L; a[2] = 6.0L;
	a[3] = 3.0L; a[4] = 5.0L; a[5] = 7.0L;
	a[6] = 4.0L; a[7] = 9.0L; a[8] = 2.0L;
}

#endif /* !MC_MAGIC3X3_H */

/* EOF */
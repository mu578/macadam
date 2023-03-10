//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_magic4x4.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_MAGIC4X4_H
#define MC_MAGIC4X4_H

#pragma mark - mc_magic4x4 -

MC_TARGET_FUNC void mc_magic4x4f(float a[16])
{
	a[0]  = 16.0f; a[1]  = 2.0f;  a[2]  = 3.0f;  a[3]  = 13.0f;
	a[4]  = 5.0f;  a[5]  = 11.0f; a[6]  = 10.0f; a[7]  = 8.0f;
	a[8]  = 9.0f;  a[9]  = 7.0f;  a[10] = 6.0f;  a[11] = 12.0f;
	a[12] = 4.0f;  a[13] = 14.0f; a[14] = 15.0f; a[15] = 1.0f;
}

MC_TARGET_FUNC void mc_magic4x4(double a[16])
{
	a[0]  = 16.0; a[1]  = 2.0;  a[2]  = 3.0;  a[3]  = 13.0;
	a[4]  = 5.0;  a[5]  = 11.0; a[6]  = 10.0; a[7]  = 8.0;
	a[8]  = 9.0;  a[9]  = 7.0;  a[10] = 6.0;  a[11] = 12.0;
	a[12] = 4.0;  a[13] = 14.0; a[14] = 15.0; a[15] = 1.0;
}

MC_TARGET_FUNC void mc_magic4x4l(long double a[16])
{
	a[0]  = 16.0L; a[1]  = 2.0L;  a[2]  = 3.0L;  a[3]  = 13.0L;
	a[4]  = 5.0L;  a[5]  = 11.0L; a[6]  = 10.0L; a[7]  = 8.0L;
	a[8]  = 9.0L;  a[9]  = 7.0L;  a[10] = 6.0L;  a[11] = 12.0L;
	a[12] = 4.0L;  a[13] = 14.0L; a[14] = 15.0L; a[15] = 1.0L;
}

#endif /* !MC_MAGIC4X4_H */

/* EOF */
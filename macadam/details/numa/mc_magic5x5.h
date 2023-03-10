//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_magic5x5.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_MAGIC5X5_H
#define MC_MAGIC5X5_H

#pragma mark - mc_magic5x5 -

MC_TARGET_FUNC void mc_magic5x5f(float a[25])
{
	a[0]  = 17.0f; a[1]  = 24.0f; a[2]  = 1.0f;  a[3]  = 8.0f;  a[4]  = 15.0f;
	a[5]  = 23.0f; a[6]  = 5.0f;  a[7]  = 7.0f;  a[8]  = 14.0f; a[9]  = 16.0f;
	a[10] = 4.0f;  a[11] = 6.0f;  a[12] = 13.0f; a[13] = 20.0f; a[14] = 22.0f;
	a[15] = 10.0f; a[16] = 12.0f; a[17] = 19.0f; a[18] = 21.0f; a[19] = 3.0f;
	a[20] = 11.0f; a[21] = 18.0f; a[22] = 25.0f; a[23] = 2.0f;  a[24] = 9.0f;
}

MC_TARGET_FUNC void mc_magic5x5(double a[25])
{
	a[0]  = 17.0; a[1]  = 24.0; a[2]  = 1.0;  a[3]  = 8.0;  a[4]  = 15.0;
	a[5]  = 23.0; a[6]  = 5.0;  a[7]  = 7.0;  a[8]  = 14.0; a[9]  = 16.0;
	a[10] = 4.0;  a[11] = 6.0;  a[12] = 13.0; a[13] = 20.0; a[14] = 22.0;
	a[15] = 10.0; a[16] = 12.0; a[17] = 19.0; a[18] = 21.0; a[19] = 3.0;
	a[20] = 11.0; a[21] = 18.0; a[22] = 25.0; a[23] = 2.0;  a[24] = 9.0;
}

MC_TARGET_FUNC void mc_magic5x5l(long double a[25])
{
	a[0]  = 17.0L; a[1]  = 24.0L; a[2]  = 1.0L;  a[3]  = 8.0L;  a[4]  = 15.0L;
	a[5]  = 23.0L; a[6]  = 5.0L;  a[7]  = 7.0L;  a[8]  = 14.0L; a[9]  = 16.0L;
	a[10] = 4.0L;  a[11] = 6.0L;  a[12] = 13.0L; a[13] = 20.0L; a[14] = 22.0L;
	a[15] = 10.0L; a[16] = 12.0L; a[17] = 19.0L; a[18] = 21.0L; a[19] = 3.0L;
	a[20] = 11.0L; a[21] = 18.0L; a[22] = 25.0L; a[23] = 2.0L;  a[24] = 9.0L;
}

#endif /* !MC_MAGIC5X5_H */

/* EOF */
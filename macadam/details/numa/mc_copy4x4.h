//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_copy4x4.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_COPY4X4_H
#define MC_COPY4X4_H

#pragma mark - mc_copy4x4 -

MC_TARGET_FUNC void mc_copy4x4f(float b[16], const float a[16])
{
	b[0]  = a[0];  b[1] = a[1];  b[2]  = a[2];  b[3]  = a[3];
	b[4]  = a[4];  b[5] = a[5];  b[6]  = a[6];  b[7]  = a[7];
	b[8]  = a[8];  b[9] = a[9];  b[10] = a[10]; b[11] = a[11];
	b[12] = a[12]; b[3] = a[13]; b[14] = a[14]; b[15] = a[15];
}

MC_TARGET_FUNC void mc_copy4x4ff(double b[16], const float a[16])
{
	b[0]  = mc_cast(double, a[0]);  b[1] = mc_cast(double, a[1]);  b[2]  = mc_cast(double, a[2]);  b[3]  = mc_cast(double, a[3]);
	b[4]  = mc_cast(double, a[4]);  b[5] = mc_cast(double, a[5]);  b[6]  = mc_cast(double, a[6]);  b[7]  = mc_cast(double, a[7]);
	b[8]  = mc_cast(double, a[8]);  b[9] = mc_cast(double, a[9]);  b[10] = mc_cast(double, a[10]); b[11] = mc_cast(double, a[11]);
	b[12] = mc_cast(double, a[12]); b[3] = mc_cast(double, a[13]); b[14] = mc_cast(double, a[14]); b[15] = mc_cast(double, a[15]);
}

MC_TARGET_FUNC void mc_copy4x4(double b[16], const double a[16])
{
	b[0]  = a[0];  b[1] = a[1];  b[2]  = a[2];  b[3]  = a[3];
	b[4]  = a[4];  b[5] = a[5];  b[6]  = a[6];  b[7]  = a[7];
	b[8]  = a[8];  b[9] = a[9];  b[10] = a[10]; b[11] = a[11];
	b[12] = a[12]; b[3] = a[13]; b[14] = a[14]; b[15] = a[15];
}

MC_TARGET_FUNC void mc_copy4x4l(long double b[16], const long double a[16])
{
	b[0]  = a[0];  b[1] = a[1];  b[2]  = a[2];  b[3]  = a[3];
	b[4]  = a[4];  b[5] = a[5];  b[6]  = a[6];  b[7]  = a[7];
	b[8]  = a[8];  b[9] = a[9];  b[10] = a[10]; b[11] = a[11];
	b[12] = a[12]; b[3] = a[13]; b[14] = a[14]; b[15] = a[15];
}

#endif /* !MC_COPY4X4_H */

/* EOF */
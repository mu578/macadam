//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_copy3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_COPY3X3_H
#define MC_COPY3X3_H

#pragma mark - mc_copy3x3 -

MC_TARGET_FUNC void mc_copy3x3f(float b[9], const float a[9])
{
	b[0] = a[0]; b[1] = a[1]; b[2] = a[2];
	b[3] = a[3]; b[4] = a[4]; b[5] = a[5];
	b[6] = a[6]; b[7] = a[7]; b[8] = a[8];
}

MC_TARGET_FUNC void mc_copy3x3ff(double b[9], const float a[9])
{
	b[0] = mc_cast(double, a[0]); b[1] = mc_cast(double, a[1]); b[2] = mc_cast(double, a[2]);
	b[3] = mc_cast(double, a[3]); b[4] = mc_cast(double, a[4]); b[5] = mc_cast(double, a[5]);
	b[6] = mc_cast(double, a[6]); b[7] = mc_cast(double, a[7]); b[8] = mc_cast(double, a[8]);
}

MC_TARGET_FUNC void mc_copy3x3(double b[9], const double a[9])
{
	b[0] = a[0]; b[1] = a[1]; b[2] = a[2];
	b[3] = a[3]; b[4] = a[4]; b[5] = a[5];
	b[6] = a[6]; b[7] = a[7]; b[8] = a[8];
}

MC_TARGET_FUNC void mc_copy3x3l(long double b[9], const long double a[9])
{
	b[0] = a[0]; b[1] = a[1]; b[2] = a[2];
	b[3] = a[3]; b[4] = a[4]; b[5] = a[5];
	b[6] = a[6]; b[7] = a[7]; b[8] = a[8];
}

#endif /* !MC_COPY3X3_H */

/* EOF */
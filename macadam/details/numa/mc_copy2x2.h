//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_copy2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_COPY2X2_H
#define MC_COPY2X2_H

#pragma mark - mc_copy2x2 -

MC_TARGET_FUNC void mc_copy2x2f(float b[4], const float a[4])
{
	b[0] = a[0]; b[1] = a[1];
	b[2] = a[2]; b[3] = a[3];
}

MC_TARGET_FUNC void mc_copy2x2ff(double b[4], const float a[4])
{
	b[0] = mc_cast(double, a[0]); b[1] = mc_cast(double, a[1]);
	b[2] = mc_cast(double, a[2]); b[3] = mc_cast(double, a[3]);
}

MC_TARGET_FUNC void mc_copy2x2(double b[4], const double a[4])
{
	b[0] = a[0]; b[1] = a[1];
	b[2] = a[2]; b[3] = a[3];
}

MC_TARGET_FUNC void mc_copy2x2l(long double b[4], const long double a[4])
{
	b[0] = a[0]; b[1] = a[1];
	b[2] = a[2]; b[3] = a[3];
}

#endif /* !MC_COPY2X2_H */

/* EOF */
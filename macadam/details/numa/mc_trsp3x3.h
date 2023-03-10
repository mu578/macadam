//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_trsp3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_trsi3x3.h>

#ifndef MC_TRSP3X3_H
#define MC_TRSP3X3_H

#pragma mark - mc_trsp3x3 -

MC_TARGET_FUNC void mc_trsp3x3f(float at[9], const float a[9])
{
//!# Returning transpose of A.
	if (a != at) {
		at[0] = a[0]; at[1] = a[3]; at[2] = a[6];
		at[3] = a[1]; at[4] = a[4]; at[5] = a[7];
		at[6] = a[2]; at[7] = a[5]; at[8] = a[8];
	} else {
		mc_trsi3x3f(at);
	}
}

MC_TARGET_FUNC void mc_trsp3x3(double at[9], const double a[9])
{
//!# Returning transpose of A.
	if (a != at) {
		at[0] = a[0]; at[1] = a[3]; at[2] = a[6];
		at[3] = a[1]; at[4] = a[4]; at[5] = a[7];
		at[6] = a[2]; at[7] = a[5]; at[8] = a[8];
	} else {
		mc_trsi3x3(at);
	}
}

MC_TARGET_FUNC void mc_trsp3x3l(long double at[9], const long double a[9])
{
//!# Returning transpose of A.
	if (a != at) {
		at[0] = a[0]; at[1] = a[3]; at[2] = a[6];
		at[3] = a[1]; at[4] = a[4]; at[5] = a[7];
		at[6] = a[2]; at[7] = a[5]; at[8] = a[8];
	} else {
		mc_trsi3x3l(at);
	}
}

#endif /* !MC_TRSP3X3_H */

/* EOF */
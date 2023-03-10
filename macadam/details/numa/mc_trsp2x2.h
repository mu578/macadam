//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_trsp2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_trsi2x2.h>

#ifndef MC_TRSP2X2_H
#define MC_TRSP2X2_H

#pragma mark - mc_trsp2x2 -

MC_TARGET_FUNC void mc_trsp2x2f(float at[4], const float a[4])
{
//!# Returning transpose of A.
	if (a != at) {
		at[0] = a[0]; at[1] = a[2];
		at[2] = a[1]; at[3] = a[3];
	} else {
		mc_trsi2x2f(at);
	}
}

MC_TARGET_FUNC void mc_trsp2x2(double at[4], const double a[4])
{
//!# Returning transpose of A.
	if (a != at) {
		at[0] = a[0]; at[1] = a[2];
		at[2] = a[1]; at[3] = a[3];
	} else {
		mc_trsi2x2(at);
	}
}

MC_TARGET_FUNC void mc_trsp2x2l(long double at[4], const long double a[4])
{
//!# Returning transpose of A.
	if (a != at) {
		at[0] = a[0]; at[1] = a[2];
		at[2] = a[1]; at[3] = a[3];
	} else {
		mc_trsi2x2l(at);
	}
}

#endif /* !MC_TRSP2X2_H */

/* EOF */
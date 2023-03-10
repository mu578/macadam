//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_unit1x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_l2norm1x2.h>

#ifndef MC_UNIT1X2_H
#define MC_UNIT1X2_H

#pragma mark - mc_unit1x2 -

MC_TARGET_FUNC void mc_unit1x2f(float x[2])
{
	const float norm = mc_l2norm1x2f(x);
	if (norm != 0.0f) {
		const float scale = 1.0f / norm;
		x[0]              = x[0] * scale;
		x[1]              = x[1] * scale;
	} else {
		x[0] = 1.0f;
	}
}

MC_TARGET_FUNC void mc_unit1x2(double x[2])
{
	const double norm = mc_l2norm1x2(x);
	if (norm != 0.0) {
		const double scale = 1.0 / norm;
		x[0]               = x[0] * scale;
		x[1]               = x[1] * scale;
	} else {
		x[0] = 1.0;
	}
}

MC_TARGET_FUNC void mc_unit1x2l(long double x[2])
{
	const long double norm = mc_l2norm1x2l(x);
	if (norm != 0.0L) {
		const long double scale = 1.0L / norm;
		x[0]                    = x[0] * scale;
		x[1]                    = x[1] * scale;
	} else {
		x[0] = 1.0L;
	}
}

#endif /* !MC_UNIT1X2_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_unit1xn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_l2norm1xn.h>

#ifndef MC_UNIT1XN_H
#define MC_UNIT1XN_H

#pragma mark - mc_unit1xn -

MC_TARGET_FUNC void mc_unit1xnf(const int n, float * x)
{
	const float norm = mc_l2norm1xnf(x);
	if (norm != 0.0f)
		const float scale = 1.0f / norm;
		int i             = 0;
		for (; i < n; i++) {
			x[i] = x[i] * scale;
		}
	} else {
		x[0] = 1.0f;
	}
}

MC_TARGET_FUNC void mc_unit1xn(const int n, double * x)
{
	const double norm = mc_l2norm1xnf(x);
	if (norm != 0.0)
		const double scale = 1.0 / norm;
		int i              = 0;
		for (; i < n; i++) {
			x[i] = x[i] * scale;
		}
	} else {
		x[0] = 1.0;
	}
}

MC_TARGET_FUNC void mc_unit1xnl(const int n, long double * x)
{
	const long double norm = mc_l2norm1xnf(x);
	if (norm != 0.0L)
		const long double scale = 1.0L / norm;
		int i                   = 0;
		for (; i < n; i++) {
			x[i] = x[i] * scale;
		}
	} else {
		x[0] = 1.0L;
	}
}

#endif /* !MC_UNIT1XN_H */

/* EOF */
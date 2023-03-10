//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lchoose.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_gammaln.h>

#ifndef MC_LCHOOSE_H
#define MC_LCHOOSE_H

#pragma mark - mc_lchoose -

MC_TARGET_FUNC float mc_lchoosef(const float x, const float y)
{
	if (0 <= y && y <= x) {
		return mc_gammalnf(x + 1.0f) - mc_gammalnf(y + 1.0f) - mc_gammalnf(x - y + 1.0f);
	}
	return MCLIMITS_MAXF;
}

MC_TARGET_FUNC double mc_lchoose(const double x, const double y)
{
	if (0 <= y && y <= x) {
		return mc_gammaln(x + 1.0) - mc_gammaln(y + 1.0) - mc_gammaln(x - y + 1.0);
	}
	return MCLIMITS_MAX;
}

MC_TARGET_FUNC long double mc_lchoosel(const long double x, const long double y)
{
	if (0 <= y && y <= x) {
		return mc_gammalnl(x + 1.0L) - mc_gammalnl(y + 1.0L) - mc_gammalnl(x - y + 1.0L);
	}
	return MCLIMITS_MAXL;
}

#endif /* !MC_LCHOOSE_H */

/* EOF */
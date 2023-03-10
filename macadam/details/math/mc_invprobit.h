//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_invprobit.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_erfc.h>

#ifndef MC_INVPROBIT_H
#define MC_INVPROBIT_H

#pragma mark - mc_invprobit -

MC_TARGET_FUNC float mc_invprobitf(const float x)
{
	return 0.5f * mc_erfcf(-x * MCK_KF(MCK_1_SQRT2));
}

MC_TARGET_FUNC double mc_invprobit(const double x)
{
	return 0.5 * mc_erfc(-x * MCK_K(MCK_1_SQRT2));
}

MC_TARGET_FUNC long double mc_invprobitl(const long double x)
{
	return 0.5L * mc_erfcl(-x * MCK_KL(MCK_1_SQRT2));
}

#endif /* !MC_INVPROBIT_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_riemann_zeta_n.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_RIEMANN_ZETA_N_H
#define MC_RIEMANN_ZETA_N_H

#pragma mark - mc_riemann_zeta_n -

MC_TARGET_FUNC float mc_riemann_zeta_nf(const float x)
{
	return x + MCK_INFN;
}

MC_TARGET_FUNC double mc_riemann_zeta_n(const double x)
{
	return x + MCK_INFN;
}

MC_TARGET_FUNC long double mc_riemann_zeta_nl(const long double x)
{
	return x + MCK_INFN;
}

#endif /* !MC_RIEMANN_ZETA_N_H */

/* EOF */
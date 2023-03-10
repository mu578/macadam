//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_riemann_zeta.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_riemann_zeta_n.h>
#include <macadam/details/math/mc_riemann_zeta_p.h>

#ifndef MC_RIEMANN_ZETA_H
#define MC_RIEMANN_ZETA_H

#pragma mark - mc_riemann_zeta -

MC_TARGET_FUNC float mc_riemann_zetaf(const float x)
{
	if (x >= 0.0f) {
		return mc_riemann_zeta_pf(x);
	}
	return mc_riemann_zeta_nf(x);
}

MC_TARGET_FUNC double mc_riemann_zeta(const double x)
{
	if (x >= 0.0) {
		return mc_riemann_zeta_p(x);
	}
	return mc_riemann_zeta_n(x);
}

MC_TARGET_FUNC long double mc_riemann_zetal(const long double x)
{
	if (x >= 0.0L) {
		return mc_riemann_zeta_pl(x);
	}
	return mc_riemann_zeta_nl(x);
}

#endif /* !MC_RIEMANN_ZETA_H */

/* EOF */
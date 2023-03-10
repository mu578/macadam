//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_log1pe.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_exp.h>
#include <macadam/details/math/mc_log1p.h>

#ifndef MC_LOG1PE_H
#define MC_LOG1PE_H

#pragma mark - mc_log1pe -

MC_TARGET_FUNC float mc_log1pef(const float x)
{
	return mc_log1pf(mc_expf(x));
}

MC_TARGET_FUNC double mc_log1pe(const double x)
{
	return mc_log1p(mc_exp(x));
}

MC_TARGET_FUNC long double mc_log1pel(const long double x)
{
	return mc_log1pl(mc_expl(x));
}

#endif /* !MC_LOG1PE_H */

/* EOF */
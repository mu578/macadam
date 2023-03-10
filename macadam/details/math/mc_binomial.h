//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_binomial.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_choose.h>

#ifndef MC_BINOMIAL_H
#define MC_BINOMIAL_H

MC_TARGET_FUNC float mc_binomialf(const unsigned int n, const unsigned int k)
{
	return mc_choosef(n, k);
}

MC_TARGET_FUNC double mc_binomial(const unsigned int n, const unsigned int k)
{
	return mc_choose(n, k);
}

MC_TARGET_FUNC long double mc_binomiall(const unsigned int n, const unsigned int k)
{
	return mc_choosel(n, k);
}

#endif /* !MC_BINOMIAL_H */

/* EOF */
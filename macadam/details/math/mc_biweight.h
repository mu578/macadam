//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_biweight.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_bisquare.h>

#ifndef MC_BIWEIGHT_H
#define MC_BIWEIGHT_H

#pragma mark - mc_biweight -

MC_TARGET_FUNC float mc_biweightf(const float r)
{
	const float c = 4.685f;
	const float s = 1.0f;
	return mc_bisquaref(r, c, s);
}

MC_TARGET_FUNC double mc_biweight(const double r)
{
	const double c = 4.685;
	const double s = 1.0;
	return mc_bisquare(r, c, s);
}

MC_TARGET_FUNC long double mc_biweightl(const long double r)
{
	const long double c = 4.685L;
	const long double s = 1.0L;
	return mc_bisquarel(r, c, s);
}

#endif /* !MC_BIWEIGHT_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_l2norm1x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_sqrt.h>
#include <macadam/details/numa/mc_ssqr1x3.h>

#ifndef MC_NORM1X3_H
#define MC_NORM1X3_H

#pragma mark - mc_l2norm1x3 -

MC_TARGET_PROC float mc_l2norm1x3f(const float x[3])
{
#	if MC_TARGET_EMBEDDED
	return mc_sqrtf(mc_raise2f(x[0]) + mc_raise2f(x[1]) + mc_raise2f(x[2]));
#	else
	float sumsq = 0.0f, scale;
	mc_ssqr1x3f(x, &sumsq, &scale);
	sumsq = scale * mc_sqrtf(sumsq);
	return sumsq;
#	endif
}

MC_TARGET_PROC double mc_l2norm1x3ff(const float x[3])
{
#	if MC_TARGET_EMBEDDED
	return mc_sqrt(mc_raise2(mc_cast(double, x[0])) + mc_raise2(mc_cast(double, x[1])) + mc_cast(double, mc_raise2(x[2])));
#	else
	double sumsq = 0.0, scale;
	mc_ssqr1x3ff(x, &sumsq, &scale);
	sumsq = scale * mc_sqrt(sumsq);
	return sumsq;
#	endif
}

MC_TARGET_PROC double mc_l2norm1x3(const double x[3])
{
#	if MC_TARGET_EMBEDDED
	return mc_sqrt(mc_raise2(x[0]) + mc_raise2(x[1]) + mc_raise2(x[2]));
#	else
	double sumsq = 0.0, scale;
	mc_ssqr1x3(x, &sumsq, &scale);
	sumsq = scale * mc_sqrt(sumsq);
	return sumsq;
#	endif
}

MC_TARGET_PROC long double mc_l2norm1x3l(const long double x[3])
{
#	if MC_TARGET_EMBEDDED
	return mc_sqrtl(mc_raise2l(x[0]) + mc_raise2l(x[1]) + mc_raise2l(x[2]));
#	else
	long double sumsq = 0.0L, scale;
	mc_ssqr1x3l(x, &sumsq, &scale);
	sumsq = scale * mc_sqrtl(sumsq);
	return sumsq;
#	endif
}

#endif /* !MC_NORM1X3_H */

/* EOF */
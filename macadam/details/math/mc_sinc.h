//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sinc.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_sin.h>

#ifndef MC_SINC_H
#define MC_SINC_H

#pragma mark - mc_sinc -

MC_TARGET_FUNC float mc_sincf(const float x)
{
//!# \note: f(0)=1, i.e removable singularity.
	const float pix = MCK_KF(MCK_PI) * x;
	return x == 0 ? 1.0f : mc_sinf(pix) / pix;
}

MC_TARGET_FUNC double mc_sinc(const double x)
{
//!# \note: f(0)=1, i.e removable singularity.
	const double pix = MCK_K(MCK_PI) * x;
	return x == 0 ? 1.0 : mc_sin(pix) / pix;
}

MC_TARGET_FUNC long double mc_sincl(const long double x)
{
//!# \note: f(0)=1, i.e removable singularity.
#	if (MC_TARGET_C99 || MC_TARGET_CPP17) && defined(M_PIl)
	const long double pix = M_PIl * x;
#	else
	const long double pix = MCK_KL(MCK_PI) * x;
#	endif
	return x == 0 ? 1.0L : mc_sinl(pix) / pix;
}

#pragma mark - mc_unnsinc -

MC_TARGET_PROC float mc_unnsincf(const float x)
{
//!# \note: f(0)=1, i.e removable singularity.
	return x == 0 ? 1.0f : mc_sinf(x) / x;
}

MC_TARGET_PROC double mc_unnsinc(const double x)
{
//!# \note: f(0)=1, i.e removable singularity.
	return x == 0 ? 1.0 : mc_sin(x) / x;
}

MC_TARGET_PROC long double mc_unnsincl(const long double x)
{
//!# \note: f(0)=1, i.e removable singularity.
	return x == 0 ? 1.0L : mc_sinl(x) / x;
}

#endif /* !MC_SINC_H */

/* EOF */
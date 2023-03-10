//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sincos.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_SINCOS_H
#define MC_SINCOS_H

#pragma mark - mc_sincos -

MC_TARGET_FUNC void mc_sincosf(const float x, float * sinp, float * cosp)
{
#	if MC_TARGET_EMBEDDED
#	if MC_TARGET_CPP98
	*sinp = ::sinf(x);
	*cosp = ::sqrtf(1.0f - (*sinp) * (*sinp));
#	else
	*sinp = sinf(x);
	*cosp = sqrtf(1.0f - (*sinp) * (*sinp));
#	endif
#	elif MC_TARGET_APPLEXM
#	if MC_TARGET_CPP98
		::__sincosf(x, sinp, cosp);
#	else
		__sincosf(x, sinp, cosp);
#	endif
#	elif (defined(__unix__) || defined(__bsdi__)) && !(defined(__linux__) || defined(__gnu_linux__))
#	if MC_TARGET_CPP98
	::sincosf(x, sinp, cosp);
#	else
	sincosf(x, sinp, cosp);
#	endif
#	else
#	if MC_TARGET_CPP98
	*sinp = ::sinf(x);
	*cosp = ::cosf(x);
#	else
	*sinp = sinf(x);
	*cosp = cosf(x);
#	endif
#	endif
}

MC_TARGET_FUNC void mc_sincos(const double x, double * sinp, double * cosp)
{
#	if MC_TARGET_EMBEDDED
#	if MC_TARGET_CPP98
	*sinp = ::sin(x);
	*cosp = ::sqrt(1.0 - (*sinp) * (*sinp));
#	else
	*sinp = sin(x);
	*cosp = sqrt(1.0 - (*sinp) * (*sinp));
#	endif
#	elif MC_TARGET_APPLEXM
#	if MC_TARGET_CPP98
	::__sincos(x, sinp, cosp);
#	else
	__sincos(x, sinp, cosp);
#	endif
#	elif (defined(__unix__) || defined(__bsdi__)) && !defined(__linux__)
#	if MC_TARGET_CPP98
	::sincos(x, sinp, cosp);
#	else
	sincos(x, sinp, cosp);
#	endif
#	else
#	if MC_TARGET_CPP98
	*sinp = ::sin(x);
	*cosp = ::cos(x);
#	else
	*sinp = sin(x);
	*cosp = cos(x);
#	endif
#	endif
}

MC_TARGET_FUNC void mc_sincosl(const long double x, long double * sinp, long double * cosp)
{
#	if MC_TARGET_EMBEDDED
#	if MC_TARGET_CPP98
	*sinp = ::sinl(x);
	*cosp = ::sqrtl(1.0L - (*sinp) * (*sinp));
#	else
	*sinp = sinl(x);
	*cosp = sqrtl(1.0L - (*sinp) * (*sinp));
#	endif
#	elif MC_TARGET_APPLEXM
	const double xx = mc_cast(double, x);
	double ss, cc;
#	if MC_TARGET_CPP98
		::__sincos(xx, &ss, &cc);
#	else
		__sincos(xx, &ss, &cc);
#	endif
	*sinp = mc_cast(long double, ss);
	*cosp = mc_cast(long double, cc);
#	elif (defined(__unix__) || defined(__bsdi__)) && !defined(__linux__)
#	if MC_TARGET_CPP98
	::sincosl(x, sinp, cosp);
#	else
	sincosl(x, sinp, cosp);
#	endif
#	else
#	if MC_TARGET_CPP98
	*sinp = ::sinl(x);
	*cosp = ::cosl(x);
#	else
	*sinp = sinl(x);
	*cosp = cosl(x);
#	endif
#	endif
}

#endif /* !MC_SINCOS_H */

/* EOF */
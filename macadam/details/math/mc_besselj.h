//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_besselj.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_BESSELJ_H
#define MC_BESSELJ_H

#pragma mark - mc_besselj0 -

MC_TARGET_FUNC float mc_besselj0f(const float x)
{
#	if  defined(__unix__)      \
	||  defined(__linux__)     \
	||  defined(__gnu_linux__) \
	||  defined(__bsdi__)      \
	||  defined(__FreeBSD__)   \
	||  defined(__NetBSD__)    \
	||  defined(__OpenBSD__)   \
	||  defined(__DragonFly__) \
	|| !defined(__APPLE__)
#	if MC_TARGET_CPP98
	return ::j0f(x);
#	else
	return j0f(x);
#	endif
#	elif MC_TARGET_MSVC_CPP
	return mc_cast(float, ::_j0(mc_cast(double, x)));
#	else
#	if MC_TARGET_CPP98
	return mc_cast(float, ::j0(mc_cast(double, x)));
#	else
	return mc_cast(float, j0(mc_cast(double, x)));
#	endif
#	endif
}

MC_TARGET_FUNC double mc_besselj0(const double x)
{
#	if MC_TARGET_CPP98
	return ::j0(x);
#	else
	return j0(x);
#	endif
}

MC_TARGET_FUNC long double mc_besselj0l(const long double x)
{
#	if  defined(__unix__)      \
	||  defined(__linux__)     \
	||  defined(__gnu_linux__) \
	||  defined(__bsdi__)      \
	||  defined(__FreeBSD__)   \
	||  defined(__NetBSD__)    \
	||  defined(__OpenBSD__)   \
	||  defined(__DragonFly__) \
	|| !defined(__APPLE__)
#	if MC_TARGET_CPP98
	return ::j0l(x);
#	else
	return j0l(x);
#	endif
#	elif MC_TARGET_MSVC_CPP
	return mc_cast(long double, ::_j0(mc_cast(double, x)));
#	else
#	if MC_TARGET_CPP98
	return mc_cast(long double, ::j0(mc_cast(double, x)));
#	else
	return mc_cast(long double, j0(mc_cast(double, x)));
#	endif
#	endif
}

#pragma mark - mc_besselj1 -

MC_TARGET_FUNC float mc_besselj1f(const float x)
{
#	if  defined(__unix__)      \
	||  defined(__linux__)     \
	||  defined(__gnu_linux__) \
	||  defined(__bsdi__)      \
	||  defined(__FreeBSD__)   \
	||  defined(__NetBSD__)    \
	||  defined(__OpenBSD__)   \
	||  defined(__DragonFly__) \
	|| !defined(__APPLE__)
#	if MC_TARGET_CPP98
	return ::j1f(x);
#	else
	return j1f(x);
#	endif
#	elif MC_TARGET_MSVC_CPP
	return mc_cast(float, ::_j1(mc_cast(double, x)));
#	else
#	if MC_TARGET_CPP98
	return mc_cast(float, ::j1(mc_cast(double, x)));
#	else
	return mc_cast(float, j1(mc_cast(double, x)));
#	endif
#	endif
}

MC_TARGET_FUNC double mc_besselj1(const double x)
{
#	if MC_TARGET_CPP98
	return ::j1(x);
#	else
	return j1(x);
#	endif
}

MC_TARGET_FUNC long double mc_besselj1l(const long double x)
{
#	if  defined(__unix__)      \
	||  defined(__linux__)     \
	||  defined(__gnu_linux__) \
	||  defined(__bsdi__)      \
	||  defined(__FreeBSD__)   \
	||  defined(__NetBSD__)    \
	||  defined(__OpenBSD__)   \
	||  defined(__DragonFly__) \
	|| !defined(__APPLE__)
#	if MC_TARGET_CPP98
	return ::j1l(x);
#	else
	return j1l(x);
#	endif
#	elif MC_TARGET_MSVC_CPP
	return mc_cast(long double, ::_j1(mc_cast(double, x)));
#	else
#	if MC_TARGET_CPP98
	return mc_cast(long double, ::j1(mc_cast(double, x)));
#	else
	return mc_cast(long double, j1(mc_cast(double, x)));
#	endif
#	endif
}

#pragma mark - mc_besseljn -

MC_TARGET_FUNC float mc_besseljnf(const int n, float x)
{
#	if  defined(__unix__)      \
	||  defined(__linux__)     \
	||  defined(__gnu_linux__) \
	||  defined(__bsdi__)      \
	||  defined(__FreeBSD__)   \
	||  defined(__NetBSD__)    \
	||  defined(__OpenBSD__)   \
	||  defined(__DragonFly__) \
	|| !defined(__APPLE__)
#	if MC_TARGET_CPP98
	return ::jnf(n, x);
#	else
	return jnf(n, x);
#	endif
#	elif MC_TARGET_MSVC_CPP
	return mc_cast(float, ::_jn(n, mc_cast(double, x)));
#	else
#	if MC_TARGET_CPP98
	return mc_cast(float, ::jn(n, mc_cast(double, x)));
#	else
	return mc_cast(float, jn(n, mc_cast(double, x)));
#	endif
#	endif
}

MC_TARGET_FUNC double mc_besseljn(const int n, double x)
{
#	if MC_TARGET_CPP98
	return ::jn(n, x);
#	else
	return jn(n, x);
#	endif
}

MC_TARGET_FUNC long double mc_besseljnl(const int n, long double x)
{
#	if  defined(__unix__)      \
	||  defined(__linux__)     \
	||  defined(__gnu_linux__) \
	||  defined(__bsdi__)      \
	||  defined(__FreeBSD__)   \
	||  defined(__NetBSD__)    \
	||  defined(__OpenBSD__)   \
	||  defined(__DragonFly__) \
	|| !defined(__APPLE__)
#	if MC_TARGET_CPP98
	return ::jnl(n, x);
#	else
	return jnl(n, x);
#	endif
#	elif MC_TARGET_MSVC_CPP
	return mc_cast(long double, ::_jn(n, mc_cast(double, x)));
#	else
#	if MC_TARGET_CPP98
	return mc_cast(long double, ::jn(n, mc_cast(double, x)));
#	else
	return mc_cast(long double, jn(n, mc_cast(double, x)));
#	endif
#	endif
}

#endif /* !MC_BESSELJ_H */

/* EOF */
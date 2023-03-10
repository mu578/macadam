//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_absmag.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/mcconsts.h>
#include <macadam/mclimits.h>

#ifndef MC_ABSMAG_H
#define MC_ABSMAG_H

#pragma mark - mc_babs -

MC_TARGET_PROC signed char mc_babs(const signed char x)
{
	return x < 0 ? -x : x;
}

#pragma mark - mc_sabs -

MC_TARGET_PROC short mc_sabs(const short x)
{
	return x < 0 ? -x : x;
}

#pragma mark - mc_iabs -

MC_TARGET_PROC int mc_iabs(const int x)
{
#	if MC_TARGET_CPP98
	return ::abs(x);
#	else
	return abs(x);
#	endif
}

#pragma mark - mc_iabs8 -

MC_TARGET_PROC int8_t mc_iabs8(const int8_t x)
{
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	return x == INT8_MIN ? INT8_MAX : (x >= 0 ? x : -x);
#	else
	return x >= 0 ? x : -x;
#	endif
}

#pragma mark - mc_iabs16 -

MC_TARGET_PROC int16_t mc_iabs16(const int16_t x)
{
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	return x == INT16_MIN ? INT16_MAX : (x >= 0 ? x : -x);
#	else
	return x >= 0 ? x : -x;
#	endif
}

#pragma mark - mc_iabs32 -

MC_TARGET_PROC int32_t mc_iabs32(const int32_t x)
{
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	return x == INT32_MIN ? INT32_MAX : (x >= 0 ? x : -x);
#	else
	return x >= 0 ? x : -x;
#	endif
}

#pragma mark - mc_iabs64 -

MC_TARGET_PROC int64_t mc_iabs64(const int64_t x)
{
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	return x == INT64_MIN ? INT64_MIN : (x >= 0 ? x : -x);
#	else
	return x >= 0 ? x : -x;
#	endif
}

#pragma mark - mc_labs -

MC_TARGET_PROC long mc_labs(const long x)
{
#	if MC_TARGET_CPP98
	return ::labs(x);
#	else
	return labs(x);
#	endif
}

#pragma mark - mc_llabs -

#	if MC_TARGET_CPP11
MC_TARGET_PROC long long mc_llabs(const long long x)
{
	return ::llabs(x);
}
#	elif MC_TARGET_C99
MC_TARGET_PROC long long mc_llabs(const long long x)
{
	return llabs(x);
}
#	else
#	define mc_llabs mc_labs
#	endif

#	if defined(__clang__)
#	pragma clang diagnostic push
#	pragma clang diagnostic ignored "-Wsign-compare"
#	define mc_absmag(x) ((x) == 0 ? 0 : ((x) < 0 ? (-(x)) : (x)))
#	pragma clang diagnostic pop
#	elif defined(__GNUC__)
#	pragma GCC diagnostic push
#	pragma GCC diagnostic ignored "-Wsign-compare"
#	define mc_absmag(x) ((x) == 0 ? 0 : ((x) < 0 ? (-(x)) : (x)))
#	pragma GCC diagnostic pop
#	elif defined(_MSC_VER)
#	pragma warning(push)
#	pragma warning(disable:4018)
#	define mc_absmag(x) ((x) == 0 ? 0 : ((x) < 0 ? (-(x)) : (x)))
#	pragma warning(pop)
#	else
#	define mc_absmag(x) ((x) == 0 ? 0 : ((x) < 0 ? (-(x)) : (x)))
#	endif

#endif /* !MC_ABSMAG_H */

/* EOF */
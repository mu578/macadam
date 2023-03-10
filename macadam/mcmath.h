//
// # -*- coding: utf-8, tab-width: 3 -*-

// mcmath.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>
#include <macadam/details/mc_math.h>

#ifndef MCMATH_H
#define MCMATH_H

#	if MC_TARGET_CPP98
#	if MC_TARGET_CPP11
#		define mcmath_float_t       ::std::float_t
#		define mcmath_double_t      ::std::double_t
#		define mcmath_fpclassify    ::std::fpclassify
#		define mcmath_isnormal      ::std::isnormal
#		define mcmath_isfinite      ::std::isfinite
#		define mcmath_isinf         ::std::isinf
#		define mcmath_isnan         ::std::isnan
#		define mcmath_signbit       ::std::signbit
#	else
#		define mcmath_float_t       ::float_t
#		define mcmath_double_t      ::double_t
#		define mcmath_fpclassify    ::fpclassify
#		define mcmath_isnormal      ::isnormal
#		define mcmath_isfinite      ::isfinite
#		define mcmath_isinf         ::isinf
#		define mcmath_isnan         ::isnan
#		define mcmath_signbit       ::signbit
#	endif
#	elif MC_TARGET_C99
#		define mcmath_float_t       float_t
#		define mcmath_double_t      double_t
#		define mcmath_fpclassify    fpclassify
#		define mcmath_isnormal      isnormal
#		define mcmath_isfinite      isfinite
#		define mcmath_isinf         isinf
#		define mcmath_isnan         isnan
#		define mcmath_signbit       signbit
#	else
#		define mcmath_float_t       float
#		define mcmath_double_t      double
#		define mcmath_fpclassify(x) 0
#		define mcmath_isnormal(x)   0
#		define mcmath_isfinite(x)   0
#		define mcmath_isinf(x)      0
#		define mcmath_isnan(x)      0
#		define mcmath_signbit(x)    0
#	endif

#	if MC_TARGET_CPP98 && MC_TARGET_ALLOW_CPP_CMATH
#		define mcmath_acos(x)          ::std::acos(x)
#		define mcmath_asin(x)          ::std::asin(x)
#		define mcmath_atan(x)          ::std::atan(x)
#		define mcmath_atan2(y, x)      ::std::atan2(y, x)
#		define mcmath_cos(x)           ::std::cos(x)
#		define mcmath_sin(x)           ::std::sin(x)
#		define mcmath_tan(x)           ::std::tan(x)
#		define mcmath_acosh(x)         ::std::acosh(x)
#		define mcmath_asinh(x)         ::std::asinh(x)
#		define mcmath_atanh(x)         ::std::atanh(x)
#		define mcmath_cosh(x)          ::std::cosh(x)
#		define mcmath_sinh(x)          ::std::sinh(x)
#		define mcmath_tanh(x)          ::std::tanh(x)
#		define mcmath_exp(x)           ::std::exp(x)
#		define mcmath_exp2(x)          ::std::exp2(x)
#		define mcmath_expm1(x)         ::std::expm1(x)
#		define mcmath_log(x)           ::std::log(x)
#		define mcmath_log10(x)         ::std::log10(x)
#		define mcmath_log2(x)          ::std::log2(x)
#	ifndef __APPLE__
#		define mcmath_log1p(x)         ::std::log1p(x)
#	endif
#		define mcmath_logb(x)          ::std::logb(x)
#		define mcmath_ldexp(x, y)      ::std::ldexp(x, y)
#		define mcmath_frexp(x, y)      ::std::frexp(x, y)
#		define mcmath_ilogb(x)         ::std::ilogb(x)
#		define mcmath_scalbn(x, y)     ::std::scalbn(x, y)
#		define mcmath_scalbln(x, y)    ::std::scalbln(x, y)
#	if MC_TARGET_CPP11
#		define mcmath_fabs(x)          ::std::fabs(x)
#		define mcmath_cbrt(x)          ::std::cbrt(x)
#		define mcmath_hypot(x, y)      ::std::hypot(x, y)
#		define mcmath_pow(x, y)        ::std::pow(x, y)
#		define mcmath_sqrt(x)          ::std::sqrt(x)
#	endif
#		define mcmath_erf(x)           ::std::erf(x)
#		define mcmath_erfc(x)          ::std::erfc(x)
#		define mcmath_lgamma(x)        ::std::lgamma(x)
#		define mcmath_tgamma(x)        ::std::tgamma(x)
#		define mcmath_ceil(x)          ::std::ceil(x)
#		define mcmath_floor(x)         ::std::floor(x)
#		define mcmath_nearbyint(x)     ::std::nearbyint(x)
#		define mcmath_rint(x)          ::std::rint(x)
#		define mcmath_lrint(x)         ::std::lrint(x)
#	if MC_TARGET_CPP11
#		define mcmath_llrint(x)        ::std::llrint(x)
#	else
#		define mcmath_llrint(x)        ::std::lrint(x)
#	endif
#		define mcmath_round(x)         ::std::round(x)
#		define mcmath_lround(x)        ::std::lround(x)
#	if MC_TARGET_CPP11
#		define mcmath_llround(x)       ::std::llround(x)
#	else
#		define mcmath_llround(x)       ::std::lround(x)
#	endif
#		define mcmath_trunc(x)         ::std::trunc(x)
#		define mcmath_fmod(x, y)       ::std::fmod(x, y)
#		define mcmath_remainder(x, y)  ::std::remainder(x, y)
#		define mcmath_remquo(x, y, q)  ::std::remquo(x, y, q)
#		define mcmath_copysign(x, y)   ::std::copysign(x, y)
#		define mcmath_nextafter(x, y)  ::std::nextafter(x, y)
#		define mcmath_nexttoward(x, y) ::std::nexttoward(x, y)
#		define mcmath_fdim(x, y)       ::std::fdim(x, y)
#	if MC_TARGET_CPP11
#		define mcmath_fmax(x, y)       ::std::fmax(x, y)
#		define mcmath_fmin(x, y)       ::std::fmin(x, y)
#	endif
#		define mcmath_fma(x, y, z)     ::std::fma(x, y, z)
#	elif MC_TARGET_C89 && MC_TARGET_HAVE_TGMATH
#		define mcmath_acos(x)          acos(x)
#		define mcmath_asin(x)          asin(x)
#		define mcmath_atan(x)          atan(x)
#		define mcmath_atan2(y, x)      atan2(y, x)
#		define mcmath_cos(x)           cos(x)
#		define mcmath_sin(x)           sin(x)
#		define mcmath_tan(x)           tan(x)
#		define mcmath_acosh(x)         acosh(x)
#		define mcmath_asinh(x)         asinh(x)
#		define mcmath_atanh(x)         atanh(x)
#		define mcmath_cosh(x)          cosh(x)
#		define mcmath_sinh(x)          sinh(x)
#		define mcmath_tanh(x)          tanh(x)
#		define mcmath_exp(x)           exp(x)
#		define mcmath_exp2(x)          exp2(x)
#		define mcmath_expm1(x)         expm1(x)
#		define mcmath_log(x)           log(x)
#		define mcmath_log10(x)         log10(x)
#		define mcmath_log2(x)          log2(x)
#	ifndef __APPLE__
#		define mcmath_log1p(x)         log1p(x)
#	endif
#		define mcmath_logb(x)          logb(x)
#		define mcmath_ldexp(x, y)      ldexp(x, y)
#		define mcmath_frexp(x, y)      frexp(x, y)
#		define mcmath_ilogb(x)         ilogb(x)
#		define mcmath_scalbn(x, y)     scalbn(x, y)
#		define mcmath_scalbln(x, y)    scalbln(x, y)
#		define mcmath_fabs(x)          fabs(x)
#		define mcmath_cbrt(x)          cbrt(x)
#		define mcmath_hypot(x, y)      hypot(x, y)
#		define mcmath_pow(x, y)        pow(x, y)
#		define mcmath_sqrt(x)          sqrt(x)
#		define mcmath_erf(x)           erf(x)
#		define mcmath_erfc(x)          erfc(x)
#		define mcmath_lgamma(x)        lgamma(x)
#		define mcmath_tgamma(x)        tgamma(x)
#		define mcmath_ceil(x)          ceil(x)
#		define mcmath_floor(x)         floor(x)
#		define mcmath_nearbyint(x)     nearbyint(x)
#		define mcmath_rint(x)          rint(x)
#		define mcmath_lrint(x)         lrint(x)
#		define mcmath_llrint(x)        llrint(x)
#		define mcmath_round(x)         round(x)
#		define mcmath_lround(x)        lround(x)
#	if MC_TARGET_C11
#		define mcmath_llround(x)       llround(x)
#	else
#		define mcmath_llround(x)       lround(x)
#	endif
#		define mcmath_trunc(x)         trunc(x)
#		define mcmath_fmod(x, y)       fmod(x, y)
#		define mcmath_remainder(x, y)  remainder(x, y)
#		define mcmath_remquo(x, y, q)  remquo(x, y, q)
#		define mcmath_copysign(x, y)   copysign(x, y)
#		define mcmath_nextafter(x, y)  nextafter(x, y)
#		define mcmath_nexttoward(x, y) nexttoward(x, y)
#		define mcmath_fdim(x, y)       fdim(x, y)
#		define mcmath_fmax(x, y)       fmax(x, y)
#		define mcmath_fmin(x, y)       fmin(x, y)
#		define mcmath_fma(x, y, z)     fma(x, y, z)
#	endif


#pragma mark - mcmath_ispow2 -

#	ifndef mcmath_ispow2
#	if MC_TARGET_CPP20 && MC_TARGET_HAVE_ISPOW2FN
#	define mcmath_ispow2(x) ::std::ispow2(x)
#	else
#	define mcmath_ispow2(x) ((x) & ((x) - 1)) == 0 && ((x) != 0)
#	endif
#	endif

#pragma mark - mcmath_nlz -

#	ifndef mcmath_nlz
#	if MC_TARGET_C99 || MC_TARGET_CPP11
#	define mcmath_nlz(x) ((x) ? MC_TARGET_CLZLL(mc_cast(unsigned long long, x)) : CHAR_BIT * sizeof(unsigned long long))
#	else
#	define mcmath_nlz(x) ((x) ? MC_TARGET_CLZL(mc_cast(unsigned long, x))       : CHAR_BIT * sizeof(unsigned long)     )
#	endif
#	endif

#pragma mark - mcmath_log2floor -

#	ifndef mcmath_log2floor
#	define mcmath_log2floor(x) ((CHAR_BIT * sizeof(x)) - 1) - mcmath_nlz(x))
#	endif

#pragma mark - mcmath_log2ceil -

#	ifndef mcmath_log2ceil
#	define mcmath_log2ceil(x) ((x) ? (CHAR_BIT * sizeof(x)) - mcmath_nlz(((x) - 1)) : -1)
#	endif

#pragma mark - mcmath_floor2 -

#	ifndef mcmath_floor2
#	if MC_TARGET_CPP20 && MC_TARGET_HAVE_FLOOR2FN
#	define mcmath_floor2(x) ::std::floor2(x)
#	else
#	if MC_TARGET_C99 || MC_TARGET_CPP11
#	define mcmath_floor2(x) ((x) ? 1ULL << mcmath_log2floor(x) : 0ULL)
#	else
#	define mcmath_floor2(x) ((x) ? 1UL  << mcmath_log2floor(x) : 0UL )
#	endif
#	endif
#	endif

#pragma mark - mcmath_ceil2 -

#	ifndef mcmath_ceil2
#	if MC_TARGET_CPP20 && MC_TARGET_HAVE_CEIL2FN
#	define mcmath_ceil2(x) ::std::ceil2(x)
#	else
#	if MC_TARGET_C99 || MC_TARGET_CPP11
#	define mcmath_ceil2(x) ((x) > 1 ? 1ULL << mcmath_log2ceil(x) : 0ULL)
#	else
#	define mcmath_ceil2(x) ((x) > 1 ? 1UL  << mcmath_log2ceil(x) : 0UL )
#	endif
#	endif
#	endif

#pragma mark - mcmath_acos -

#	ifndef mcmath_acos
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_acos              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_acos<float>       (const float& x)       { return mc_acosf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_acos<double>      (const double& x)      { return mc_acos  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_acos<long double> (const long double& x) { return mc_acosl (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_acos (const float x)       { return mc_acosf (x); }
MC_TARGET_ALIAS double      mcmath_acos (const double x)      { return mc_acos  (x); }
MC_TARGET_ALIAS long double mcmath_acos (const long double x) { return mc_acosl (x); }
#	elif MC_TARGET_C11
#	define mcmath_acos(x) _Generic(x \
	, float       : mc_acosf         \
	, double      : mc_acos          \
	, long double : mc_acosl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_acos(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_acosf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_acos  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_acosl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_acos(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_acosf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_acos  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_acosl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_asin -

#	ifndef mcmath_asin
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_asin              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_asin<float>       (const float& x)       { return mc_asinf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_asin<double>      (const double& x)      { return mc_asin  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_asin<long double> (const long double& x) { return mc_asinl (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_asin (const float x)       { return mc_asinf (x); }
MC_TARGET_ALIAS double      mcmath_asin (const double x)      { return mc_asin  (x); }
MC_TARGET_ALIAS long double mcmath_asin (const long double x) { return mc_asinl (x); }
#	elif MC_TARGET_C11
#	define mcmath_asin(x) _Generic(x \
	, float       : mc_asinf         \
	, double      : mc_asin          \
	, long double : mc_asinl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_asin(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_asinf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_asin  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_asinl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_asin(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_asinf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_asin  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_asinl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_atan -

#	ifndef mcmath_atan
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_atan              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_atan<float>       (const float& x)       { return mc_atanf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_atan<double>      (const double& x)      { return mc_atan  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_atan<long double> (const long double& x) { return mc_atanl (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_atan (const float x)       { return mc_atanf (x); }
MC_TARGET_ALIAS double      mcmath_atan (const double x)      { return mc_atan  (x); }
MC_TARGET_ALIAS long double mcmath_atan (const long double x) { return mc_atanl (x); }
#	elif MC_TARGET_C11
#	define mcmath_atan(x) _Generic(x \
	, float       : mc_atanf         \
	, double      : mc_atan          \
	, long double : mc_atanl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_atan(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_atanf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_atan  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_atanl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_atan(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_atanf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_atan  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_atanl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_atan2 -

#	ifndef mcmath_atan2
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_atan2              (const T& y, const T& x)                     { mc_unused(y); mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_atan2<float>       (const float& y, const float& x)             { return mc_atan2f (y, x);              }
template <>        MC_TARGET_INLINE double      mcmath_atan2<double>      (const double& y, const double& x)           { return mc_atan2  (y, x);              }
template <>        MC_TARGET_INLINE long double mcmath_atan2<long double> (const long double& y, const long double& x) { return mc_atan2l (y, x);              }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_atan2 (float y, float x)             { return mc_atan2f (y, x); }
MC_TARGET_ALIAS double      mcmath_atan2 (double y, double x)           { return mc_atan2  (y, x); }
MC_TARGET_ALIAS long double mcmath_atan2 (long double y, long double x) { return mc_atan2l (y, x); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_atan2(y, x) _Generic(x \
	, float       : mc_atan2f            \
	, double      : mc_atan2             \
	, long double : mc_atan2l            \
) (y, mc_cast_expr(MC_TARGET_TYPEOF(y), x))
#	elif MC_TARGET_C11
#	define mcmath_atan2(y, x) _Generic((y)+(x) \
	, float       : mc_atan2f                  \
	, double      : mc_atan2                   \
	, long double : mc_atan2l                  \
) ((y), (x))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_atan2(y, x) mc_cast(MC_TARGET_TYPEOF(y), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(y), float)       ? mc_atan2f (mc_cast_expr(const float, y), mc_cast_expr(const float, x))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(y), double)      ? mc_atan2  (mc_cast_expr(const double, y), mc_cast_expr(const double, x))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(y), long double) ? mc_atan2l (mc_cast_expr(const long double, y), mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_atan2(y, x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_atan2f (mc_cast_expr(const float, y), mc_cast_expr(const float, x))             \
		: sizeof(x) == sizeof(double)      ? mc_atan2  (mc_cast_expr(const double, y), mc_cast_expr(const double, x))           \
		: sizeof(x) == sizeof(long double) ? mc_atan2l (mc_cast_expr(const long double, y), mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_cos -

#	ifndef mcmath_cos
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_cos              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_cos<float>       (const float& x)       { return mc_cosf (x);     }
template <>        MC_TARGET_INLINE double      mcmath_cos<double>      (const double& x)      { return mc_cos  (x);     }
template <>        MC_TARGET_INLINE long double mcmath_cos<long double> (const long double& x) { return mc_cosl (x);     }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_cos (const float x)       { return mc_cosf (x); }
MC_TARGET_ALIAS double      mcmath_cos (const double x)      { return mc_cos  (x); }
MC_TARGET_ALIAS long double mcmath_cos (const long double x) { return mc_cosl (x); }
#	elif MC_TARGET_C11
#	define mcmath_cos(x) _Generic(x \
	, float       : mc_cosf         \
	, double      : mc_cos          \
	, long double : mc_cosl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_cos(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_cosf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_cos  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_cosl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_cos(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_cosf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_cos  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_cosl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_cospi -

#	ifndef mcmath_cospi
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_cospi              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_cospi<float>       (const float& x)       { return mc_cospif (x);   }
template <>        MC_TARGET_INLINE double      mcmath_cospi<double>      (const double& x)      { return mc_cospi  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_cospi<long double> (const long double& x) { return mc_cospil (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_cospi (const float x)       { return mc_cospif (x); }
MC_TARGET_ALIAS double      mcmath_cospi (const double x)      { return mc_cospi  (x); }
MC_TARGET_ALIAS long double mcmath_cospi (const long double x) { return mc_cospil (x); }
#	elif MC_TARGET_C11
#	define mcmath_cospi(x) _Generic(x \
	, float       : mc_cospif         \
	, double      : mc_cospi          \
	, long double : mc_cospil         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_cospi(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_cospif (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_cospi  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_cospil (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_cospi(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_cospif (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_cospi  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_cospil (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_sin -

#	ifndef mcmath_sin
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_sin              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_sin<float>       (const float& x)       { return mc_sinf (x);     }
template <>        MC_TARGET_INLINE double      mcmath_sin<double>      (const double& x)      { return mc_sin  (x);     }
template <>        MC_TARGET_INLINE long double mcmath_sin<long double> (const long double& x) { return mc_sinl (x);     }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_sin (const float x)       { return mc_sinf (x); }
MC_TARGET_ALIAS double      mcmath_sin (const double x)      { return mc_sin  (x); }
MC_TARGET_ALIAS long double mcmath_sin (const long double x) { return mc_sinl (x); }
#	elif MC_TARGET_C11
#	define mcmath_sin(x) _Generic(x \
	, float       : mc_sinf         \
	, double      : mc_sin          \
	, long double : mc_sinl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_sin(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_sinf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_sin  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_sinl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_sin(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_sinf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_sin  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_sinl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_sinpi -

#	ifndef mcmath_sinpi
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_sinpi              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_sinpi<float>       (const float& x)       { return mc_sinpif (x);   }
template <>        MC_TARGET_INLINE double      mcmath_sinpi<double>      (const double& x)      { return mc_sinpi  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_sinpi<long double> (const long double& x) { return mc_sinpil (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_sinpi (const float x)       { return mc_sinpif (x); }
MC_TARGET_ALIAS double      mcmath_sinpi (const double x)      { return mc_sinpi  (x); }
MC_TARGET_ALIAS long double mcmath_sinpi (const long double x) { return mc_sinpil (x); }
#	elif MC_TARGET_C11
#	define mcmath_sinpi(x) _Generic(x \
	, float       : mc_sinpif         \
	, double      : mc_sinpi          \
	, long double : mc_sinpil         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_sinpi(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_sinpif (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_sinpi  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_sinpil (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_sinpi(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_sinpif (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_sinpi  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_sinpil (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_tan -

#	ifndef mcmath_tan
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_tan              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_tan<float>       (const float& x)       { return mc_tanf (x);     }
template <>        MC_TARGET_INLINE double      mcmath_tan<double>      (const double& x)      { return mc_tan  (x);     }
template <>        MC_TARGET_INLINE long double mcmath_tan<long double> (const long double& x) { return mc_tanl (x);     }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_tan (const float x)       { return mc_tanf (x); }
MC_TARGET_ALIAS double      mcmath_tan (const double x)      { return mc_tan  (x); }
MC_TARGET_ALIAS long double mcmath_tan (const long double x) { return mc_tanl (x); }
#	elif MC_TARGET_C11
#	define mcmath_tan(x) _Generic(x \
	, float       : mc_tanf         \
	, double      : mc_tan          \
	, long double : mc_tanl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_tan(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_tanf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_tan  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_tanl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_tan(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_tanf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_tan  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_tanl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_acosh -

#	ifndef mcmath_acosh
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_acosh              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_acosh<float>       (const float& x)       { return mc_acoshf (x);   }
template <>        MC_TARGET_INLINE double      mcmath_acosh<double>      (const double& x)      { return mc_acosh  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_acosh<long double> (const long double& x) { return mc_acoshl (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_acosh (const float x)       { return mc_acoshf (x); }
MC_TARGET_ALIAS double      mcmath_acosh (const double x)      { return mc_acosh  (x); }
MC_TARGET_ALIAS long double mcmath_acosh (const long double x) { return mc_acoshl (x); }
#	elif MC_TARGET_C11
#	define mcmath_acosh(x) _Generic(x \
	, float       : mc_acoshf         \
	, double      : mc_acosh          \
	, long double : mc_acoshl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_acosh(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_acoshf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_acosh  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_acoshl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_acosh(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_acoshf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_acosh  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_acoshl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_asinh -

#	ifndef mcmath_asinh
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_asinh              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_asinh<float>       (const float& x)       { return mc_asinhf (x);   }
template <>        MC_TARGET_INLINE double      mcmath_asinh<double>      (const double& x)      { return mc_asinh  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_asinh<long double> (const long double& x) { return mc_asinhl (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_asinh (const float x)       { return mc_asinhf (x); }
MC_TARGET_ALIAS double      mcmath_asinh (const double x)      { return mc_asinh  (x); }
MC_TARGET_ALIAS long double mcmath_asinh (const long double x) { return mc_asinhl (x); }
#	elif MC_TARGET_C11
#	define mcmath_asinh(x) _Generic(x \
	, float       : mc_asinhf         \
	, double      : mc_asinh          \
	, long double : mc_asinhl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_asinh(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_asinhf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_asinh  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_asinhl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_asinh(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_asinhf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_asinh  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_asinhl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_atanh -

#	ifndef mcmath_atanh
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_atanh              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_atanh<float>       (const float& x)       { return mc_atanhf (x);   }
template <>        MC_TARGET_INLINE double      mcmath_atanh<double>      (const double& x)      { return mc_atanh  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_atanh<long double> (const long double& x) { return mc_atanhl (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_atanh (const float x)       { return mc_atanhf (x); }
MC_TARGET_ALIAS double      mcmath_atanh (const double x)      { return mc_atanh  (x); }
MC_TARGET_ALIAS long double mcmath_atanh (const long double x) { return mc_atanhl (x); }
#	elif MC_TARGET_C11
#	define mcmath_atanh(x) _Generic(x \
	, float       : mc_atanhf         \
	, double      : mc_atanh          \
	, long double : mc_atanhl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_atanh(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_atanhf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_atanh  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_atanhl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_atanh(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_atanhf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_atanh  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_atanhl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_cosh -

#	ifndef mcmath_cosh
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_cosh              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_cosh<float>       (const float& x)       { return mc_coshf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_cosh<double>      (const double& x)      { return mc_cosh  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_cosh<long double> (const long double& x) { return mc_coshl (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_cosh (const float x)       { return mc_coshf (x); }
MC_TARGET_ALIAS double      mcmath_cosh (const double x)      { return mc_cosh  (x); }
MC_TARGET_ALIAS long double mcmath_cosh (const long double x) { return mc_coshl (x); }
#	elif MC_TARGET_C11
#	define mcmath_cosh(x) _Generic(x \
	, float       : mc_coshf         \
	, double      : mc_cosh          \
	, long double : mc_coshl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_cosh(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_coshf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_cosh  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_coshl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_cosh(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_coshf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_cosh  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_coshl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_sinh -

#	ifndef mcmath_sinh
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_sinh              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_sinh<float>       (const float& x)       { return mc_sinhf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_sinh<double>      (const double& x)      { return mc_sinh  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_sinh<long double> (const long double& x) { return mc_sinhl (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_sinh (const float x)       { return mc_sinhf (x); }
MC_TARGET_ALIAS double      mcmath_sinh (const double x)      { return mc_sinh  (x); }
MC_TARGET_ALIAS long double mcmath_sinh (const long double x) { return mc_sinhl (x); }
#	elif MC_TARGET_C11
#	define mcmath_sinh(x) _Generic(x \
	, float       : mc_sinhf         \
	, double      : mc_sinh          \
	, long double : mc_sinhl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_sinh(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_sinhf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_sinh  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_sinhl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_sinh(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_sinhf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_sinh  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_sinhl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_tanh -

#	ifndef mcmath_tanh
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_tanh              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_tanh<float>       (const float& x)       { return mc_tanhf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_tanh<double>      (const double& x)      { return mc_tanh  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_tanh<long double> (const long double& x) { return mc_tanhl (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_tanh (const float x)       { return mc_tanhf (x); }
MC_TARGET_ALIAS double      mcmath_tanh (const double x)      { return mc_tanh  (x); }
MC_TARGET_ALIAS long double mcmath_tanh (const long double x) { return mc_tanhl (x); }
#	elif MC_TARGET_C11
#	define mcmath_tanh(x) _Generic(x \
	, float       : mc_tanhf         \
	, double      : mc_tanh          \
	, long double : mc_tanhl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_tanh(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_tanhf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_tanh  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_tanhl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_tanh(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_tanhf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_tanh  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_tanhl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_exp -

#	ifndef mcmath_exp
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_exp              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_exp<float>       (const float& x)       { return mc_expf (x);     }
template <>        MC_TARGET_INLINE double      mcmath_exp<double>      (const double& x)      { return mc_exp  (x);     }
template <>        MC_TARGET_INLINE long double mcmath_exp<long double> (const long double& x) { return mc_expl (x);     }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_exp (const float x)       { return mc_expf (x); }
MC_TARGET_ALIAS double      mcmath_exp (const double x)      { return mc_exp  (x); }
MC_TARGET_ALIAS long double mcmath_exp (const long double x) { return mc_expl (x); }
#	elif MC_TARGET_C11
#	define mcmath_exp(x) _Generic(x \
	, float       : mc_expf         \
	, double      : mc_exp          \
	, long double : mc_expl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_exp(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_expf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_exp  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_expl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_exp(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_expf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_exp  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_expl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_exp2 -

#	ifndef mcmath_exp2
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_exp2              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_exp2<float>       (const float& x)       { return mc_exp2f (x);    }
template <>        MC_TARGET_INLINE double      mcmath_exp2<double>      (const double& x)      { return mc_exp2  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_exp2<long double> (const long double& x) { return mc_exp2l (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_exp2 (const float x)       { return mc_exp2f (x); }
MC_TARGET_ALIAS double      mcmath_exp2 (const double x)      { return mc_exp2  (x); }
MC_TARGET_ALIAS long double mcmath_exp2 (const long double x) { return mc_exp2l (x); }
#	elif MC_TARGET_C11
#	define mcmath_exp2(x) _Generic(x \
	, float       : mc_exp2f         \
	, double      : mc_exp2          \
	, long double : mc_exp2l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_exp2(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_exp2f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_exp2  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_exp2l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_exp2(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_exp2f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_exp2  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_exp2l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_exp2m1 -

#	ifndef mcmath_exp2m1
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_exp2m1              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_exp2m1<float>       (const float& x)       { return mc_exp2m1f (x);  }
template <>        MC_TARGET_INLINE double      mcmath_exp2m1<double>      (const double& x)      { return mc_exp2m1  (x);  }
template <>        MC_TARGET_INLINE long double mcmath_exp2m1<long double> (const long double& x) { return mc_exp2m1l (x);  }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_exp2m1 (const float x)       { return mc_exp2m1f (x); }
MC_TARGET_ALIAS double      mcmath_exp2m1 (const double x)      { return mc_exp2m1  (x); }
MC_TARGET_ALIAS long double mcmath_exp2m1 (const long double x) { return mc_exp2m1l (x); }
#	elif MC_TARGET_C11
#	define mcmath_exp2m1(x) _Generic(x \
	, float       : mc_exp2m1f         \
	, double      : mc_exp2m1          \
	, long double : mc_exp2m1l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_exp2m1(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_exp2m1f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_exp2m1  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_exp2m1l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_exp2m1(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_exp2m1f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_exp2m1  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_exp2m1l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_exp10 -

#	ifndef mcmath_exp10
#	if MC_TARGET_CPP98

template <class T> MC_TARGET_INLINE T           mcmath_exp10              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_exp10<float>       (const float& x)       { return mc_exp10f (x);   }
template <>        MC_TARGET_INLINE double      mcmath_exp10<double>      (const double& x)      { return mc_exp10  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_exp10<long double> (const long double& x) { return mc_exp10l (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_exp10 (const float x)       { return mc_exp10f (x); }
MC_TARGET_ALIAS double      mcmath_exp10 (const double x)      { return mc_exp10  (x); }
MC_TARGET_ALIAS long double mcmath_exp10 (const long double x) { return mc_exp10l (x); }
#	elif MC_TARGET_C11
#	define mcmath_exp10(x) _Generic(x \
	, float       : mc_exp10f         \
	, double      : mc_exp10          \
	, long double : mc_exp10l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_exp10(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_exp10f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_exp10  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_exp10l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_exp10(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_exp10f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_exp10  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_exp10l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_exp10m1 -

#	ifndef mcmath_exp10m1
#	if MC_TARGET_CPP98

template <class T> MC_TARGET_INLINE T           mcmath_exp10m1              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_exp10m1<float>       (const float& x)       { return mc_exp10m1f (x); }
template <>        MC_TARGET_INLINE double      mcmath_exp10m1<double>      (const double& x)      { return mc_exp10m1  (x); }
template <>        MC_TARGET_INLINE long double mcmath_exp10m1<long double> (const long double& x) { return mc_exp10m1l (x); }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_exp10m1 (const float x)       { return mc_exp10m1f (x); }
MC_TARGET_ALIAS double      mcmath_exp10m1 (const double x)      { return mc_exp10m1  (x); }
MC_TARGET_ALIAS long double mcmath_exp10m1 (const long double x) { return mc_exp10m1l (x); }
#	elif MC_TARGET_C11
#	define mcmath_exp10m1(x) _Generic(x \
	, float       : mc_exp10m1f         \
	, double      : mc_exp10m1          \
	, long double : mc_exp10m1l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_exp10m1(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_exp10m1f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_exp10m1  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_exp10m1l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_exp10m1(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_exp10m1f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_exp10m1  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_exp10m1l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_expit -

#	ifndef mcmath_expit
#	if MC_TARGET_CPP98

template <class T> MC_TARGET_INLINE T           mcmath_expit              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_expit<float>       (const float& x)       { return mc_expitf (x);   }
template <>        MC_TARGET_INLINE double      mcmath_expit<double>      (const double& x)      { return mc_expit  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_expit<long double> (const long double& x) { return mc_expitl (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_expit (const float x)       { return mc_expitf (x); }
MC_TARGET_ALIAS double      mcmath_expit (const double x)      { return mc_expit  (x); }
MC_TARGET_ALIAS long double mcmath_expit (const long double x) { return mc_expitl (x); }
#	elif MC_TARGET_C11
#	define mcmath_expit(x) _Generic(x \
	, float       : mc_expitf         \
	, double      : mc_expit          \
	, long double : mc_expitl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_expit(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_expitf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_expit  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_expitl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_expit(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_expitf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_expit  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_expitl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_expm1 -

#	ifndef mcmath_expm1
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_expm1              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_expm1<float>       (const float& x)       { return mc_expm1f (x);   }
template <>        MC_TARGET_INLINE double      mcmath_expm1<double>      (const double& x)      { return mc_expm1  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_expm1<long double> (const long double& x) { return mc_expm1l (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_expm1 (const float x)       { return mc_expm1f (x); }
MC_TARGET_ALIAS double      mcmath_expm1 (const double x)      { return mc_expm1  (x); }
MC_TARGET_ALIAS long double mcmath_expm1 (const long double x) { return mc_expm1l (x); }
#	elif MC_TARGET_C11
#	define mcmath_expm1(x) _Generic(x \
	, float       : mc_expm1f         \
	, double      : mc_expm1          \
	, long double : mc_expm1l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_expm1(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_expm1f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_expm1  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_expm1l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_expm1(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_expm1f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_expm1  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_expm1l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_log -

#	ifndef mcmath_log
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_log              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_log<float>       (const float& x)       { return mc_logf (x);     }
template <>        MC_TARGET_INLINE double      mcmath_log<double>      (const double& x)      { return mc_log  (x);     }
template <>        MC_TARGET_INLINE long double mcmath_log<long double> (const long double& x) { return mc_logl (x);     }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_log (const float x)       { return mc_logf (x); }
MC_TARGET_ALIAS double      mcmath_log (const double x)      { return mc_log  (x); }
MC_TARGET_ALIAS long double mcmath_log (const long double x) { return mc_logl (x); }
#	elif MC_TARGET_C11
#	define mcmath_log(x) _Generic(x \
	, float       : mc_logf         \
	, double      : mc_log          \
	, long double : mc_logl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_log(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_logf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_log  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_logl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_log(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_logf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_log  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_logl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_log10 -

#	ifndef mcmath_log10
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_log10              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_log10<float>       (const float& x)       { return mc_log10f (x);   }
template <>        MC_TARGET_INLINE double      mcmath_log10<double>      (const double& x)      { return mc_log10  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_log10<long double> (const long double& x) { return mc_log10l (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_log10 (const float x)       { return mc_log10f (x); }
MC_TARGET_ALIAS double      mcmath_log10 (const double x)      { return mc_log10  (x); }
MC_TARGET_ALIAS long double mcmath_log10 (const long double x) { return mc_log10l (x); }
#	elif MC_TARGET_C11
#	define mcmath_log10(x) _Generic(x \
	, float       : mc_log10f         \
	, double      : mc_log10          \
	, long double : mc_log10l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_log10(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_log10f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_log10  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_log10l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_log10(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_log10f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_log10  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_log10l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_log2 -

#	ifndef mcmath_log2
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_log2              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_log2<float>       (const float& x)       { return mc_log2f (x);    }
template <>        MC_TARGET_INLINE double      mcmath_log2<double>      (const double& x)      { return mc_log2  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_log2<long double> (const long double& x) { return mc_log2l (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_log2 (const float x)       { return mc_log2f (x); }
MC_TARGET_ALIAS double      mcmath_log2 (const double x)      { return mc_log2  (x); }
MC_TARGET_ALIAS long double mcmath_log2 (const long double x) { return mc_log2l (x); }
#	elif MC_TARGET_C11
#	define mcmath_log2(x) _Generic(x \
	, float       : mc_log2f         \
	, double      : mc_log2          \
	, long double : mc_log2l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_log2(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_log2f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_log2  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_log2l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_log2(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_log2f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_log2  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_log2l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_log1m -

#	ifndef mcmath_log1m
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_log1m              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_log1m<float>       (const float& x)       { return mc_log1mf (x);   }
template <>        MC_TARGET_INLINE double      mcmath_log1m<double>      (const double& x)      { return mc_log1m  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_log1m<long double> (const long double& x) { return mc_log1ml (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_log1m (const float x)       { return mc_log1mf (x); }
MC_TARGET_ALIAS double      mcmath_log1m (const double x)      { return mc_log1m  (x); }
MC_TARGET_ALIAS long double mcmath_log1m (const long double x) { return mc_log1ml (x); }
#	elif MC_TARGET_C11
#	define mcmath_log1m(x) _Generic(x \
	, float       : mc_log1mf         \
	, double      : mc_log1m          \
	, long double : mc_log1ml         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_log1m(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_log1mf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_log1m  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_log1ml (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_log1m(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_log1mf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_log1m  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_log1ml (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_log1me -

#	ifndef mcmath_log1me
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_log1me              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_log1me<float>       (const float& x)       { return mc_log1mef (x);  }
template <>        MC_TARGET_INLINE double      mcmath_log1me<double>      (const double& x)      { return mc_log1me  (x);  }
template <>        MC_TARGET_INLINE long double mcmath_log1me<long double> (const long double& x) { return mc_log1mel (x);  }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_log1me (const float x)       { return mc_log1mef (x); }
MC_TARGET_ALIAS double      mcmath_log1me (const double x)      { return mc_log1me  (x); }
MC_TARGET_ALIAS long double mcmath_log1me (const long double x) { return mc_log1mel (x); }
#	elif MC_TARGET_C11
#	define mcmath_log1me(x) _Generic(x \
	, float       : mc_log1mef         \
	, double      : mc_log1me          \
	, long double : mc_log1mel         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_log1me(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_log1mef (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_log1me  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_log1mel (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_log1me(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_log1mef (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_log1me  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_log1mel (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_log1p -

#	ifndef mcmath_log1p
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_log1p              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_log1p<float>       (const float& x)       { return mc_log1pf (x);   }
template <>        MC_TARGET_INLINE double      mcmath_log1p<double>      (const double& x)      { return mc_log1p  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_log1p<long double> (const long double& x) { return mc_log1pl (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_log1p (const float x)       { return mc_log1pf (x); }
MC_TARGET_ALIAS double      mcmath_log1p (const double x)      { return mc_log1p  (x); }
MC_TARGET_ALIAS long double mcmath_log1p (const long double x) { return mc_log1pl (x); }
#	elif MC_TARGET_C11
#	define mcmath_log1p(x) _Generic(x \
	, float       : mc_log1pf         \
	, double      : mc_log1p          \
	, long double : mc_log1pl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_log1p(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_log1pf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_log1p  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_log1pl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_log1p(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_log1pf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_log1p  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_log1pl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_log1pe -

#	ifndef mcmath_log1pe
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_log1pe              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_log1pe<float>       (const float& x)       { return mc_log1pef (x);  }
template <>        MC_TARGET_INLINE double      mcmath_log1pe<double>      (const double& x)      { return mc_log1pe  (x);  }
template <>        MC_TARGET_INLINE long double mcmath_log1pe<long double> (const long double& x) { return mc_log1pel (x);  }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_log1pe (const float x)       { return mc_log1pef (x); }
MC_TARGET_ALIAS double      mcmath_log1pe (const double x)      { return mc_log1pe  (x); }
MC_TARGET_ALIAS long double mcmath_log1pe (const long double x) { return mc_log1pel (x); }
#	elif MC_TARGET_C11
#	define mcmath_log1pe(x) _Generic(x \
	, float       : mc_log1pef         \
	, double      : mc_log1pe          \
	, long double : mc_log1pel         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_log1pe(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_log1pef (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_log1pe  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_log1pel (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_log1pe(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_log1pef (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_log1pe  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_log1pel (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_log1pmx -

#	ifndef mcmath_log1pmx
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_log1pmx              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_log1pmx<float>       (const float& x)       { return mc_log1pmxf (x); }
template <>        MC_TARGET_INLINE double      mcmath_log1pmx<double>      (const double& x)      { return mc_log1pmx  (x); }
template <>        MC_TARGET_INLINE long double mcmath_log1pmx<long double> (const long double& x) { return mc_log1pmxl (x); }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_log1pmx (const float x)       { return mc_log1pmxf (x); }
MC_TARGET_ALIAS double      mcmath_log1pmx (const double x)      { return mc_log1pmx  (x); }
MC_TARGET_ALIAS long double mcmath_log1pmx (const long double x) { return mc_log1pmxl (x); }
#	elif MC_TARGET_C11
#	define mcmath_log1pmx(x) _Generic(x \
	, float       : mc_log1pmxf         \
	, double      : mc_log1pmx          \
	, long double : mc_log1pmxl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_log1pmx(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_log1pmxf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_log1pmx  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_log1pmxl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_log1pmx(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_log1pmxf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_log1pmx  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_log1pmxl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_logp1 -

#	ifndef mcmath_logp1
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_logp1              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_logp1<float>       (const float& x)       { return mc_logp1f (x);   }
template <>        MC_TARGET_INLINE double      mcmath_logp1<double>      (const double& x)      { return mc_logp1  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_logp1<long double> (const long double& x) { return mc_logp1l (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_logp1 (const float x)       { return mc_logp1f (x); }
MC_TARGET_ALIAS double      mcmath_logp1 (const double x)      { return mc_logp1  (x); }
MC_TARGET_ALIAS long double mcmath_logp1 (const long double x) { return mc_logp1l (x); }
#	elif MC_TARGET_C11
#	define mcmath_logp1(x) _Generic(x \
	, float       : mc_logp1f         \
	, double      : mc_logp1          \
	, long double : mc_logp1l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_logp1(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_logp1f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_logp1  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_logp1l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_logp1(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_logp1f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_logp1  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_logp1l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_logb -

#	ifndef mcmath_logb
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_logb              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_logb<float>       (const float& x)       { return mc_logbf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_logb<double>      (const double& x)      { return mc_logb  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_logb<long double> (const long double& x) { return mc_logbl (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_logb (const float x)       { return mc_logbf (x); }
MC_TARGET_ALIAS double      mcmath_logb (const double x)      { return mc_logb  (x); }
MC_TARGET_ALIAS long double mcmath_logb (const long double x) { return mc_logbl (x); }
#	elif MC_TARGET_C11
#	define mcmath_logb(x) _Generic(x \
	, float       : mc_logbf         \
	, double      : mc_logb          \
	, long double : mc_logbl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_logb(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_logbf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_logb  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_logbl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_logb(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_logbf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_logb  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_logbl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_modf -

#	ifndef mcmath_modf
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_modf              (const T& x, T * y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_modf<float>       (const float& x, float * y)             { return mc_modff (x, y);               }
template <>        MC_TARGET_INLINE double      mcmath_modf<double>      (const double& x, double * y)           { return mc_modf  (x, y);               }
template <>        MC_TARGET_INLINE long double mcmath_modf<long double> (const long double& x, long double * y) { return mc_modfl (x, y);               }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_modf (const float x, float * y)             { return mc_modff (x, y); }
MC_TARGET_ALIAS double      mcmath_modf (const double x, double * y)           { return mc_modf  (x, y); }
MC_TARGET_ALIAS long double mcmath_modf (const long double x, long double * y) { return mc_modfl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_modf(x, y) _Generic(x \
	, float       : mc_modff            \
	, double      : mc_modf             \
	, long double : mc_modfl            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x) *, y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_modf(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_modff (mc_cast_expr(const float, x), mc_cast_expr(float *, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_modf  (mc_cast_expr(const double, x), mc_cast_expr(double *, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_modfl (mc_cast_expr(const long double, x), mc_cast_expr(long double *, y)) \
		: 0 \
	))
#	else
#	define mcmath_modf(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_modff (mc_cast_expr(const float, x), mc_cast_expr(float *, y))             \
		: sizeof(x) == sizeof(double)      ? mc_modf  (mc_cast_expr(const double, x), mc_cast_expr(double *, y))           \
		: sizeof(x) == sizeof(long double) ? mc_modfl (mc_cast_expr(const long double, x), mc_cast_expr(long double *, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_ldexp -

#	ifndef mcmath_ldexp
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_ldexp              (const T& x, const int n)           { mc_unused(x); mc_unused(n); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_ldexp<float>       (const float& x, const int n)       { return mc_ldexpf (x, n);              }
template <>        MC_TARGET_INLINE double      mcmath_ldexp<double>      (const double& x, const int n)      { return mc_ldexp  (x, n);              }
template <>        MC_TARGET_INLINE long double mcmath_ldexp<long double> (const long double& x, const int n) { return mc_ldexpl (x, n);              }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_ldexp (const float x, const int n)       { return mc_ldexpf (x, n); }
MC_TARGET_ALIAS double      mcmath_ldexp (const double x, const int n)      { return mc_ldexp  (x, n); }
MC_TARGET_ALIAS long double mcmath_ldexp (const long double x, const int n) { return mc_ldexpl (x, n); }
#	elif MC_TARGET_C11
#	define mcmath_ldexp(x, n) _Generic(x \
	, float       : mc_ldexpf            \
	, double      : mc_ldexp             \
	, long double : mc_ldexpl            \
) (x, mc_cast_expr(const int, n))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_ldexp(x, n) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_ldexpf (mc_cast_expr(const float, x), mc_cast_expr(const int, n))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_ldexp  (mc_cast_expr(const double, x), mc_cast_expr(const int, n))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_ldexpl (mc_cast_expr(const long double, x), mc_cast_expr(const int, n)) \
		: 0 \
	))
#	else
#	define mcmath_ldexp(x, n) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_ldexpf (mc_cast_expr(const float, x), mc_cast_expr(const int, n))       \
		: sizeof(x) == sizeof(double)      ? mc_ldexp  (mc_cast_expr(const double, x), mc_cast_expr(const int, n))      \
		: sizeof(x) == sizeof(long double) ? mc_ldexpl (mc_cast_expr(const long double, x), mc_cast_expr(const int, n)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_frexp -

#	ifndef mcmath_frexp
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_frexp              (const T& x, int * e)           { mc_unused(x); mc_unused(e); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_frexp<float>       (const float& x, int * e)       { return mc_frexpf (x, e);              }
template <>        MC_TARGET_INLINE double      mcmath_frexp<double>      (const double& x, int * e)      { return mc_frexp  (x, e);              }
template <>        MC_TARGET_INLINE long double mcmath_frexp<long double> (const long double& x, int * e) { return mc_frexpl (x, e);              }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_frexp (const float x, int * e)       { return mc_frexpf (x, e); }
MC_TARGET_ALIAS double      mcmath_frexp (const double x, int * e)      { return mc_frexp  (x, e); }
MC_TARGET_ALIAS long double mcmath_frexp (const long double x, int * e) { return mc_frexpl (x, e); }
#	elif MC_TARGET_C11
#	define mcmath_frexp(x, e) _Generic(x \
	, float       : mc_frexpf            \
	, double      : mc_frexp             \
	, long double : mc_frexpl            \
) (x, mc_cast_expr(int *, e))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_frexp(x, e) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_frexpf (mc_cast_expr(const float, x), mc_cast_expr(int *, e))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_frexp  (mc_cast_expr(const double, x), mc_cast_expr(int *, e))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_frexpl (mc_cast_expr(const long double, x), mc_cast_expr(int *, e)) \
		: 0 \
	))
#	else
#	define mcmath_frexp(x, e) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_frexpf (mc_cast_expr(const float, x), mc_cast_expr(int *, e))       \
		: sizeof(x) == sizeof(double)      ? mc_frexp  (mc_cast_expr(const double, x), mc_cast_expr(int *, e))      \
		: sizeof(x) == sizeof(long double) ? mc_frexpl (mc_cast_expr(const long double, x), mc_cast_expr(int *, e)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_ilogb -

#	ifndef mcmath_ilogb
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE int mcmath_ilogb              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE int mcmath_ilogb<float>       (const float& x)       { return mc_ilogbf (x);   }
template <>        MC_TARGET_INLINE int mcmath_ilogb<double>      (const double& x)      { return mc_ilogb  (x);   }
template <>        MC_TARGET_INLINE int mcmath_ilogb<long double> (const long double& x) { return mc_ilogbl (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS int mcmath_ilogb (const float x)       { return mc_ilogbf (x); }
MC_TARGET_ALIAS int mcmath_ilogb (const double x)      { return mc_ilogb  (x); }
MC_TARGET_ALIAS int mcmath_ilogb (const long double x) { return mc_ilogbl (x); }
#	elif MC_TARGET_C11
#	define mcmath_ilogb(x) _Generic(x \
	, float       : mc_ilogbf         \
	, double      : mc_ilogb          \
	, long double : mc_ilogbl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_ilogb(x) \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_ilogbf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_ilogb  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_ilogbl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	else
#	define mcmath_ilogb(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_ilogbf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_ilogb  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_ilogbl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_scalbn -

#	ifndef mcmath_scalbn
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_scalbn              (const T& x, const int y)           { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_scalbn<float>       (const float& x, const int y)       { return mc_scalbnf (x, y);             }
template <>        MC_TARGET_INLINE double      mcmath_scalbn<double>      (const double& x, const int y)      { return mc_scalbn  (x, y);             }
template <>        MC_TARGET_INLINE long double mcmath_scalbn<long double> (const long double& x, const int y) { return mc_scalbnl (x, y);             }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_scalbn (const float x, const int y)       { return mc_scalbnf (x, y); }
MC_TARGET_ALIAS double      mcmath_scalbn (const double x, const int y)      { return mc_scalbn  (x, y); }
MC_TARGET_ALIAS long double mcmath_scalbn (const long double x, const int y) { return mc_scalbnl (x, y); }
#	elif MC_TARGET_C11
#	define mcmath_scalbn(x, y) _Generic(x \
	, float       : mc_scalbnf            \
	, double      : mc_scalbn             \
	, long double : mc_scalbnl            \
) (x, mc_cast_expr(const int, y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_scalbn(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_scalbnf (mc_cast_expr(const float, x), mc_cast_expr(const int, y))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_scalbn  (mc_cast_expr(const double, x), mc_cast_expr(const int, y))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_scalbnl (mc_cast_expr(const long double, x), mc_cast_expr(const int, y)) \
		: 0 \
	))
#	else
#	define mcmath_scalbn(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_scalbnf (mc_cast_expr(const float, x), mc_cast_expr(const int, y))       \
		: sizeof(x) == sizeof(double)      ? mc_scalbn  (mc_cast_expr(const double, x), mc_cast_expr(const int, y))      \
		: sizeof(x) == sizeof(long double) ? mc_scalbnl (mc_cast_expr(const long double, x), mc_cast_expr(const int, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_scalbln -

#	ifndef mcmath_scalbln
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_scalbln              (const T& x, const long y)           { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_scalbln<float>       (const float& x, const long y)       { return mc_scalblnf (x, y);            }
template <>        MC_TARGET_INLINE double      mcmath_scalbln<double>      (const double& x, const long y)      { return mc_scalbln  (x, y);            }
template <>        MC_TARGET_INLINE long double mcmath_scalbln<long double> (const long double& x, const long y) { return mc_scalblnl (x, y);            }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_scalbln (const float x, const long y)       { return mc_scalblnf (x, y); }
MC_TARGET_ALIAS double      mcmath_scalbln (const double x, const long y)      { return mc_scalbln  (x, y); }
MC_TARGET_ALIAS long double mcmath_scalbln (const long double x, const long y) { return mc_scalblnl (x, y); }
#	elif MC_TARGET_C11
#	define mcmath_scalbln(x, y) _Generic(x \
	, float       : mc_scalblnf            \
	, double      : mc_scalbln             \
	, long double : mc_scalblnl            \
) (x, mc_cast_expr(const long, y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_scalbln(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_scalblnf (mc_cast_expr(const float, x), mc_cast_expr(const long, y))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_scalbln  (mc_cast_expr(const double, x), mc_cast_expr(const long, y))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_scalblnl (mc_cast_expr(const long double, x), mc_cast_expr(const long, y)) \
		: 0 \
	))
#	else
#	define mcmath_scalbln(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_scalblnf (mc_cast_expr(const float, x), mc_cast_expr(const long, y))       \
		: sizeof(x) == sizeof(double)      ? mc_scalbln  (mc_cast_expr(const double, x), mc_cast_expr(const long, y))      \
		: sizeof(x) == sizeof(long double) ? mc_scalblnl (mc_cast_expr(const long double, x), mc_cast_expr(const long, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_fabs -

#	ifndef mcmath_fabs
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_fabs              (const T& x)           { return ::std::abs(x); }
template <>        MC_TARGET_INLINE float       mcmath_fabs<float>       (const float& x)       { return mc_fabsf (x);  }
template <>        MC_TARGET_INLINE double      mcmath_fabs<double>      (const double& x)      { return mc_fabs  (x);  }
template <>        MC_TARGET_INLINE long double mcmath_fabs<long double> (const long double& x) { return mc_fabsl (x);  }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_fabs (const float x)       { return mc_fabsf (x); }
MC_TARGET_ALIAS double      mcmath_fabs (const double x)      { return mc_fabs  (x); }
MC_TARGET_ALIAS long double mcmath_fabs (const long double x) { return mc_fabsl (x); }
#	elif MC_TARGET_C11
#	define mcmath_fabs(x) _Generic(x \
	, float       : mc_fabsf         \
	, double      : mc_fabs          \
	, long double : mc_fabsl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_fabs(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_fabsf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_fabs  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_fabsl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_fabs(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_fabsf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_fabs  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_fabsl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_cbrt -

#	ifndef mcmath_cbrt
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_cbrt              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_cbrt<float>       (const float& x)       { return mc_cbrtf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_cbrt<double>      (const double& x)      { return mc_cbrt  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_cbrt<long double> (const long double& x) { return mc_cbrtl (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_cbrt (const float x)       { return mc_cbrtf (x); }
MC_TARGET_ALIAS double      mcmath_cbrt (const double x)      { return mc_cbrt  (x); }
MC_TARGET_ALIAS long double mcmath_cbrt (const long double x) { return mc_cbrtl (x); }
#	elif MC_TARGET_C11
#	define mcmath_cbrt(x) _Generic(x \
	, float       : mc_cbrtf         \
	, double      : mc_cbrt          \
	, long double : mc_cbrtl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_cbrt(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_cbrtf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_cbrt  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_cbrtl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_cbrt(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_cbrtf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_cbrt  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_cbrtl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_fhrt -

#	ifndef mcmath_fhrt
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_fhrt              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_fhrt<float>       (const float& x)       { return mc_fhrtf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_fhrt<double>      (const double& x)      { return mc_fhrt  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_fhrt<long double> (const long double& x) { return mc_fhrtl (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_fhrt (const float x)       { return mc_fhrtf (x); }
MC_TARGET_ALIAS double      mcmath_fhrt (const double x)      { return mc_fhrt  (x); }
MC_TARGET_ALIAS long double mcmath_fhrt (const long double x) { return mc_fhrtl (x); }
#	elif MC_TARGET_C11
#	define mcmath_fhrt(x) _Generic(x \
	, float       : mc_fhrtf         \
	, double      : mc_fhrt          \
	, long double : mc_fhrtl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_fhrt(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_fhrtf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_fhrt  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_fhrtl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_fhrt(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_fhrtf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_fhrt  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_fhrtl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_rootn -

#	ifndef mcmath_rootn
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_rootn              (const unsigned int n, const T& x)           { mc_unused(x); mc_unused(n); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_rootn<float>       (const unsigned int n, const float& x)       { return mc_rootnf (n, x);              }
template <>        MC_TARGET_INLINE double      mcmath_rootn<double>      (const unsigned int n, const double& x)      { return mc_rootn  (n, x);              }
template <>        MC_TARGET_INLINE long double mcmath_rootn<long double> (const unsigned int n, const long double& x) { return mc_rootnl (n, x);              }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_rootn (const unsigned int n, float x)       { return mc_rootnf (n, x); }
MC_TARGET_ALIAS double      mcmath_rootn (const unsigned int n, double x)      { return mc_rootn  (n, x); }
MC_TARGET_ALIAS long double mcmath_rootn (const unsigned int n, long double x) { return mc_rootnl (n, x); }
#	elif MC_TARGET_C11
#	define mcmath_rootn(n, x) _Generic(x \
	, float       : mc_rootnf            \
	, double      : mc_rootn             \
	, long double : mc_rootnl            \
) (mc_cast_expr(const unsigned int, n), x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_rootn(n, x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_rootnf (mc_cast_expr(const unsigned int, n), mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_rootn  (mc_cast_expr(const unsigned int, n), mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_rootnl (mc_cast_expr(const unsigned int, n), mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_rootn(n, x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_rootnf (mc_cast_expr(const unsigned int, n), mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_rootn  (mc_cast_expr(const unsigned int, n), mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_rootnl (mc_cast_expr(const unsigned int, n), mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_hypot -

#	ifndef mcmath_hypot
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_hypot              (const T& x, const T& y)                     { return mc_hypot  (mc_cast(double, x), mc_cast(double, y)); }
template <>        MC_TARGET_INLINE float       mcmath_hypot<float>       (const float& x, const float& y)             { return mc_hypotf (x, y);                                   }
template <>        MC_TARGET_INLINE double      mcmath_hypot<double>      (const double& x, const double& y)           { return mc_hypot  (x, y);                                   }
template <>        MC_TARGET_INLINE long double mcmath_hypot<long double> (const long double& x, const long double& y) { return mc_hypotl (x, y);                                   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_hypot (const float x, const float y)             { return mc_hypotf (x, y); }
MC_TARGET_ALIAS double      mcmath_hypot (const double x, const double y)           { return mc_hypot  (x, y); }
MC_TARGET_ALIAS long double mcmath_hypot (const long double x, const long double y) { return mc_hypotl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_hypot(x, y) _Generic(x \
	, float       : mc_hypotf \
	, double      : mc_hypot  \
	, long double : mc_hypotl \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_hypot(x, y) _Generic((x)+(y) \
	, float       : mc_hypotf \
	, double      : mc_hypot  \
	, long double : mc_hypotl \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_hypot(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_hypotf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_hypot  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_hypotl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_hypot(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_hypotf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_hypot  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_hypotl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#	ifndef mcmath_mc_hypot2
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_mc_hypot2              (const T& x, const T& y)                     { return mc_hypot2  (mc_cast(double, x), mc_cast(double, y)); }
template <>        MC_TARGET_INLINE float       mcmath_mc_hypot2<float>       (const float& x, const float& y)             { return mc_hypot2f (x, y);                                   }
template <>        MC_TARGET_INLINE double      mcmath_mc_hypot2<double>      (const double& x, const double& y)           { return mc_hypot2  (x, y);                                   }
template <>        MC_TARGET_INLINE long double mcmath_mc_hypot2<long double> (const long double& x, const long double& y) { return mc_hypot2l (x, y);                                   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_mc_hypot2 (const float x, const float y)             { return mc_hypot2f (x, y); }
MC_TARGET_ALIAS double      mcmath_mc_hypot2 (const double x, const double y)           { return mc_hypot2  (x, y); }
MC_TARGET_ALIAS long double mcmath_mc_hypot2 (const long double x, const long double y) { return mc_hypot2l (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_mc_hypot2(x, y) _Generic(x \
	, float       : mc_hypot2f \
	, double      : mc_hypot2  \
	, long double : mc_hypot2l \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_mc_hypot2(x, y) _Generic((x)+(y) \
	, float       : mc_hypot2f \
	, double      : mc_hypot2  \
	, long double : mc_hypot2l \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_mc_hypot2(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_hypot2f (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_hypot2  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_hypot2l (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_mc_hypot2(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_hypot2f (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_hypot2  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_hypot2l (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#	ifndef mcmath_hypot3
#	if MC_TARGET_CPP98
#	if MC_TARGET_CPP17
template <class T> MC_TARGET_INLINE T           mcmath_hypot3              (const T& x, const T& y, const T& z)                               { return ::std::hypot (mc_cast(double, x), mc_cast(double, y), mc_cast(double, z)); }
template <>        MC_TARGET_INLINE float       mcmath_hypot3<float>       (const float& x, const float& y, const float& z)                   { return ::std::hypot (x, y, z);                                                    }
template <>        MC_TARGET_INLINE double      mcmath_hypot3<double>      (const double& x, const double& y, const double& z)                { return ::std::hypot (x, y, z);                                                    }
template <>        MC_TARGET_INLINE long double mcmath_hypot3<long double> (const long double& x, const long double& y, const long double& z) { return ::std::hypot (x, y, z);                                                    }
#	else
template <class T> MC_TARGET_INLINE T           mcmath_hypot3              (const T& x, const T& y, const T& z)                               { return mc_hypot3  (mc_cast(double, x), mc_cast(double, y), mc_cast(double, z));   }
template <>        MC_TARGET_INLINE float       mcmath_hypot3<float>       (const float& x, const float& y, const float& z)                   { return mc_hypot3f (x, y, z);                                                      }
template <>        MC_TARGET_INLINE double      mcmath_hypot3<double>      (const double& x, const double& y, const double& z)                { return mc_hypot3  (x, y, z);                                                      }
template <>        MC_TARGET_INLINE long double mcmath_hypot3<long double> (const long double& x, const long double& y, const long double& z) { return mc_hypot3l (x, y, z);                                                      }
#	endif
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_hypot3 (const float x, const float y, const float z)                   { return mc_hypot3f (x, y, z); }
MC_TARGET_ALIAS double      mcmath_hypot3 (const double x, const double y, const double z)                { return mc_hypot3  (x, y, z); }
MC_TARGET_ALIAS long double mcmath_hypot3 (const long double x, const long double y, const long double z) { return mc_hypot3l (x, y, z); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_hypot3(x, y, z) _Generic(x \
	, float       : mc_hypot3f               \
	, double      : mc_hypot3                \
	, long double : mc_hypot3l               \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y), mc_cast_expr(MC_TARGET_TYPEOF(x), z))
#	elif MC_TARGET_C11
#	define mcmath_hypot3(x, y, z) _Generic((x)+(y)+(z) \
	, float       : mc_hypot3f                         \
	, double      : mc_hypot3                          \
	, long double : mc_hypot3l                         \
) ((x), (y), (z))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_hypot3(x, y, z) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_hypot3f (mc_cast_expr(const float, x), mc_cast_expr(const float, y), mc_cast_expr(const float, z))                   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_hypot3  (mc_cast_expr(const double, x), mc_cast_expr(const double, y), mc_cast_expr(const double, z))                \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_hypot3l (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y), mc_cast_expr(const long double, z)) \
		: 0 \
	))
#	else
#	define mcmath_hypot3(x, y, z) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_hypot3f (mc_cast_expr(const float, x), mc_cast_expr(const float, y), mc_cast_expr(const float, z))                   \
		: sizeof(x) == sizeof(double)      ? mc_hypot3  (mc_cast_expr(const double, x), mc_cast_expr(const double, y), mc_cast_expr(const double, z))                \
		: sizeof(x) == sizeof(long double) ? mc_hypot3l (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y), mc_cast_expr(const long double, z)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_pow -

#	ifndef mcmath_pow
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_pow              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_pow<float>       (const float& x, const float& y)             { return mc_powf (x, y);                }
template <>        MC_TARGET_INLINE double      mcmath_pow<double>      (const double& x, const double& y)           { return mc_pow  (x, y);                }
template <>        MC_TARGET_INLINE long double mcmath_pow<long double> (const long double& x, const long double& y) { return mc_powl (x, y);                }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_pow (const float x, const float y)             { return mc_powf (x, y); }
MC_TARGET_ALIAS double      mcmath_pow (const double x, const double y)           { return mc_pow  (x, y); }
MC_TARGET_ALIAS long double mcmath_pow (const long double x, const long double y) { return mc_powl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_pow(x, y) _Generic(x \
	, float       : mc_powf            \
	, double      : mc_pow             \
	, long double : mc_powl            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_pow(x, y) _Generic((x)+(y) \
	, float       : mc_powf                  \
	, double      : mc_pow                   \
	, long double : mc_powl                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_pow(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_powf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_pow  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_powl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_pow(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_powf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_pow  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_powl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_pow2 -

#	ifndef mcmath_pow2
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_pow2              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_pow2<float>       (const float& x)       { return mc_pow2f (x);    }
template <>        MC_TARGET_INLINE double      mcmath_pow2<double>      (const double& x)      { return mc_pow2  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_pow2<long double> (const long double& x) { return mc_pow2l (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_pow2 (const float x)       { return mc_pow2f (x); }
MC_TARGET_ALIAS double      mcmath_pow2 (const double x)      { return mc_pow2  (x); }
MC_TARGET_ALIAS long double mcmath_pow2 (const long double x) { return mc_pow2l (x); }
#	elif MC_TARGET_C11
#	define mcmath_pow2(x) _Generic(x \
	, float       : mc_pow2f         \
	, double      : mc_pow2          \
	, long double : mc_pow2l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_pow2(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_pow2f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_pow2  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_pow2l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_pow2(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_pow2f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_pow2  (mc_cast_expr(const double, x))       \
		: sizeof(x) == sizeof(long double) ? mc_pow2l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_sqrt -

#	ifndef mcmath_sqrt
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_sqrt              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_sqrt<float>       (const float& x)       { return mc_sqrtf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_sqrt<double>      (const double& x)      { return mc_sqrt  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_sqrt<long double> (const long double& x) { return mc_sqrtl (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_sqrt (const float x)       { return mc_sqrtf (x); }
MC_TARGET_ALIAS double      mcmath_sqrt (const double x)      { return mc_sqrt  (x); }
MC_TARGET_ALIAS long double mcmath_sqrt (const long double x) { return mc_sqrtl (x); }
#	elif MC_TARGET_C11
#	define mcmath_sqrt(x) _Generic(x \
	, float       : mc_sqrtf         \
	, double      : mc_sqrt          \
	, long double : mc_sqrtl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_sqrt(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_sqrtf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_sqrt  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_sqrtl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_sqrt(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_sqrtf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_sqrt  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_sqrtl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_sqrt1pm1 -

#	ifndef mcmath_sqrt1pm1
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_sqrt1pm1              (const T& x)           { mc_unused(x); return 0;  }
template <>        MC_TARGET_INLINE float       mcmath_sqrt1pm1<float>       (const float& x)       { return mc_sqrt1pm1f (x); }
template <>        MC_TARGET_INLINE double      mcmath_sqrt1pm1<double>      (const double& x)      { return mc_sqrt1pm1  (x); }
template <>        MC_TARGET_INLINE long double mcmath_sqrt1pm1<long double> (const long double& x) { return mc_sqrt1pm1l (x); }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_sqrt1pm1 (const float x)       { return mc_sqrt1pm1f (x); }
MC_TARGET_ALIAS double      mcmath_sqrt1pm1 (const double x)      { return mc_sqrt1pm1  (x); }
MC_TARGET_ALIAS long double mcmath_sqrt1pm1 (const long double x) { return mc_sqrt1pm1l (x); }
#	elif MC_TARGET_C11
#	define mcmath_sqrt1pm1(x) _Generic(x \
	, float       : mc_sqrt1pm1f         \
	, double      : mc_sqrt1pm1          \
	, long double : mc_sqrt1pm1l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_sqrt1pm1(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_sqrt1pm1f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_sqrt1pm1  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_sqrt1pm1l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_sqrt1pm1(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_sqrt1pm1f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_sqrt1pm1  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_sqrt1pm1l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_rsqrt -

#	ifndef mcmath_rsqrt
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_rsqrt              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_rsqrt<float>       (const float& x)       { return mc_rsqrtf (x);   }
template <>        MC_TARGET_INLINE double      mcmath_rsqrt<double>      (const double& x)      { return mc_rsqrt  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_rsqrt<long double> (const long double& x) { return mc_rsqrtl (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_rsqrt (const float x)       { return mc_rsqrtf (x); }
MC_TARGET_ALIAS double      mcmath_rsqrt (const double x)      { return mc_rsqrt  (x); }
MC_TARGET_ALIAS long double mcmath_rsqrt (const long double x) { return mc_rsqrtl (x); }
#	elif MC_TARGET_C11
#	define mcmath_rsqrt(x) _Generic(x \
	, float       : mc_rsqrtf         \
	, double      : mc_rsqrt          \
	, long double : mc_rsqrtl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_rsqrt(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_rsqrtf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_rsqrt  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_rsqrtl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_rsqrt(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_rsqrtf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_rsqrt  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_rsqrtl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_erf -

#	ifndef mcmath_erf
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_erf              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_erf<float>       (const float& x)       { return mc_erff (x);     }
template <>        MC_TARGET_INLINE double      mcmath_erf<double>      (const double& x)      { return mc_erf  (x);     }
template <>        MC_TARGET_INLINE long double mcmath_erf<long double> (const long double& x) { return mc_erfl (x);     }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_erf (const float x)       { return mc_erff (x); }
MC_TARGET_ALIAS double      mcmath_erf (const double x)      { return mc_erf  (x); }
MC_TARGET_ALIAS long double mcmath_erf (const long double x) { return mc_erfl (x); }
#	elif MC_TARGET_C11
#	define mcmath_erf(x) _Generic(x \
	, float       : mc_erff         \
	, double      : mc_erf          \
	, long double : mc_erfl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_erf(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_erff (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_erf  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_erfl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_erf(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_erff (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_erf  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_erfl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_erfc -

#	ifndef mcmath_erfc
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_erfc              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_erfc<float>       (const float& x)       { return mc_erfcf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_erfc<double>      (const double& x)      { return mc_erfc  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_erfc<long double> (const long double& x) { return mc_erfcl (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_erfc (const float x)       { return mc_erfcf (x); }
MC_TARGET_ALIAS double      mcmath_erfc (const double x)      { return mc_erfc  (x); }
MC_TARGET_ALIAS long double mcmath_erfc (const long double x) { return mc_erfcl (x); }
#	elif MC_TARGET_C11
#	define mcmath_erfc(x) _Generic(x \
	, float       : mc_erfcf         \
	, double      : mc_erfc          \
	, long double : mc_erfcl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_erfc(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_erfcf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_erfc  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_erfcl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_erfc(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_erfcf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_erfc  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_erfcl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_lgamma -

#	ifndef mcmath_lgamma
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_lgamma              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_lgamma<float>       (const float& x)       { return mc_lgammaf (x);  }
template <>        MC_TARGET_INLINE double      mcmath_lgamma<double>      (const double& x)      { return mc_lgamma  (x);  }
template <>        MC_TARGET_INLINE long double mcmath_lgamma<long double> (const long double& x) { return mc_lgammal (x);  }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_lgamma (const float x)       { return mc_lgammaf (x); }
MC_TARGET_ALIAS double      mcmath_lgamma (const double x)      { return mc_lgamma  (x); }
MC_TARGET_ALIAS long double mcmath_lgamma (const long double x) { return mc_lgammal (x); }
#	elif MC_TARGET_C11
#	define mcmath_lgamma(x) _Generic(x \
	, float       : mc_lgammaf         \
	, double      : mc_lgamma          \
	, long double : mc_lgammal         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_lgamma(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_lgammaf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_lgamma  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_lgammal (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_lgamma(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_lgammaf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_lgamma  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_lgammal (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_lgamma_r -

#	ifndef mcmath_lgamma_r
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_lgamma_r              (const T& x, int * s)           { mc_unused(x); return 0;     }
template <>        MC_TARGET_INLINE float       mcmath_lgamma_r<float>       (const float& x, int * s)       { return mc_lgammaf_r (x, s); }
template <>        MC_TARGET_INLINE double      mcmath_lgamma_r<double>      (const double& x, int * s)      { return mc_lgamma_r  (x, s); }
template <>        MC_TARGET_INLINE long double mcmath_lgamma_r<long double> (const long double& x, int * s) { return mc_lgammal_r (x, s); }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_lgamma_r (const float x, int * s)       { return mc_lgammaf_r (x, s); }
MC_TARGET_ALIAS double      mcmath_lgamma_r (const double x, int * s)      { return mc_lgamma_r  (x, s); }
MC_TARGET_ALIAS long double mcmath_lgamma_r (const long double x, int * s) { return mc_lgammal_r (x, s); }
#	elif MC_TARGET_C11
#	define mcmath_lgamma_r(x, s) _Generic(x \
	, float       : mc_lgammaf_r            \
	, double      : mc_lgamma_r             \
	, long double : mc_lgammal_r            \
) (x, mc_cast_expr(int *, s))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_lgamma_r(x, s) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_lgammaf_r (mc_cast_expr(const float, x), mc_cast_expr(int *, s))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_lgamma_r  (mc_cast_expr(const double, x), mc_cast_expr(int *, s))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_lgammal_r (mc_cast_expr(const long double, x), mc_cast_expr(int *, s)) \
		: 0 \
	))
#	else
#	define mcmath_lgamma_r(x, s) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_lgammaf_r (mc_cast_expr(const float, x), mc_cast_expr(int *, s))       \
		: sizeof(x) == sizeof(double)      ? mc_lgamma_r  (mc_cast_expr(const double, x), mc_cast_expr(int *, s))      \
		: sizeof(x) == sizeof(long double) ? mc_lgammal_r (mc_cast_expr(const long double, x), mc_cast_expr(int *, s)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_tgamma -

#	ifndef mcmath_tgamma
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_tgamma              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_tgamma<float>       (const float& x)       { return mc_tgammaf (x);  }
template <>        MC_TARGET_INLINE double      mcmath_tgamma<double>      (const double& x)      { return mc_tgamma  (x);  }
template <>        MC_TARGET_INLINE long double mcmath_tgamma<long double> (const long double& x) { return mc_tgammal (x);  }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_tgamma (const float x)       { return mc_tgammaf (x); }
MC_TARGET_ALIAS double      mcmath_tgamma (const double x)      { return mc_tgamma  (x); }
MC_TARGET_ALIAS long double mcmath_tgamma (const long double x) { return mc_tgammal (x); }
#	elif MC_TARGET_C11
#	define mcmath_tgamma(x) _Generic(x \
	, float       : mc_tgammaf         \
	, double      : mc_tgamma          \
	, long double : mc_tgammal         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_tgamma(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_tgammaf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_tgamma  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_tgammal (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_tgamma(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_tgammaf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_tgamma  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_tgammal (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_digamma -

#	ifndef mcmath_digamma
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_digamma              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_digamma<float>       (const float& x)       { return mc_digammaf (x); }
template <>        MC_TARGET_INLINE double      mcmath_digamma<double>      (const double& x)      { return mc_digamma  (x); }
template <>        MC_TARGET_INLINE long double mcmath_digamma<long double> (const long double& x) { return mc_digammal (x); }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_digamma (const float x)       { return mc_digammaf (x); }
MC_TARGET_ALIAS double      mcmath_digamma (const double x)      { return mc_digamma  (x); }
MC_TARGET_ALIAS long double mcmath_digamma (const long double x) { return mc_digammal (x); }
#	elif MC_TARGET_C11
#	define mcmath_digamma(x) _Generic(x \
	, float       : mc_digammaf         \
	, double      : mc_digamma          \
	, long double : mc_digammal         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_digamma(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_digammaf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_digamma  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_digammal (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_digamma(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_digammaf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_digamma  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_digammal (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_trigamma -

#	ifndef mcmath_trigamma
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_trigamma              (const T& x)           { mc_unused(x); return 0;  }
template <>        MC_TARGET_INLINE float       mcmath_trigamma<float>       (const float& x)       { return mc_trigammaf (x); }
template <>        MC_TARGET_INLINE double      mcmath_trigamma<double>      (const double& x)      { return mc_trigamma  (x); }
template <>        MC_TARGET_INLINE long double mcmath_trigamma<long double> (const long double& x) { return mc_trigammal (x); }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_trigamma (const float x)       { return mc_trigammaf (x); }
MC_TARGET_ALIAS double      mcmath_trigamma (const double x)      { return mc_trigamma  (x);  }
MC_TARGET_ALIAS long double mcmath_trigamma (const long double x) { return mc_trigammal (x); }
#	elif MC_TARGET_C11
#	define mcmath_trigamma(x) _Generic(x \
	, float       : mc_trigammaf         \
	, double      : mc_trigamma          \
	, long double : mc_trigammal         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_trigamma(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_trigammaf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_trigamma  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_trigammal (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_trigamma(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_trigammaf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_trigamma  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_trigammal (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#	define mcmath_gamma(x)   mcmath_tgamma(x)
#	define mcmath_gammaln(x) mcmath_lgamma(x)

#pragma mark - mcmath_ceil -

#	ifndef mcmath_ceil
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_ceil              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_ceil<float>       (const float& x)       { return mc_ceilf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_ceil<double>      (const double& x)      { return mc_ceil  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_ceil<long double> (const long double& x) { return mc_ceill (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_ceil (const float x)       { return mc_ceilf (x); }
MC_TARGET_ALIAS double      mcmath_ceil (const double x)      { return mc_ceil  (x); }
MC_TARGET_ALIAS long double mcmath_ceil (const long double x) { return mc_ceill (x); }
#	elif MC_TARGET_C11
#	define mcmath_ceil(x) _Generic(x \
	, float       : mc_ceilf         \
	, double      : mc_ceil          \
	, long double : mc_ceill         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_ceil(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_ceilf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_ceil  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_ceill (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_ceil(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_ceilf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_ceil  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_ceill (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_floor -

#	ifndef mcmath_floor
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_floor              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_floor<float>       (const float& x)       { return mc_floorf (x);   }
template <>        MC_TARGET_INLINE double      mcmath_floor<double>      (const double& x)      { return mc_floor  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_floor<long double> (const long double& x) { return mc_floorl (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_floor (const float x)       { return mc_floorf (x); }
MC_TARGET_ALIAS double      mcmath_floor (const double x)      { return mc_floor  (x); }
MC_TARGET_ALIAS long double mcmath_floor (const long double x) { return mc_floorl (x); }
#	elif MC_TARGET_C11
#	define mcmath_floor(x) _Generic(x \
	, float       : mc_floorf         \
	, double      : mc_floor          \
	, long double : mc_floorl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_floor(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_floorf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_floor  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_floorl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_floor(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_floorf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_floor  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_floorl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_nearbyint -

#	ifndef mcmath_nearbyint
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_nearbyint              (const T& x)           { mc_unused(x); return 0;   }
template <>        MC_TARGET_INLINE float       mcmath_nearbyint<float>       (const float& x)       { return mc_nearbyintf (x); }
template <>        MC_TARGET_INLINE double      mcmath_nearbyint<double>      (const double& x)      { return mc_nearbyint  (x); }
template <>        MC_TARGET_INLINE long double mcmath_nearbyint<long double> (const long double& x) { return mc_nearbyintl (x); }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_nearbyint (const float x)       { return mc_nearbyintf (x); }
MC_TARGET_ALIAS double      mcmath_nearbyint (const double x)      { return mc_nearbyint  (x); }
MC_TARGET_ALIAS long double mcmath_nearbyint (const long double x) { return mc_nearbyintl (x); }
#	elif MC_TARGET_C11
#	define mcmath_nearbyint(x) _Generic(x \
	, float       : mc_nearbyintf         \
	, double      : mc_nearbyint          \
	, long double : mc_nearbyintl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_nearbyint(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_nearbyintf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_nearbyint  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_nearbyintl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_nearbyint(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_nearbyintf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_nearbyint  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_nearbyintl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_rint -

#	ifndef mcmath_rint
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_rint              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_rint<float>       (const float& x)       { return mc_rintf (x);    }
template <>        MC_TARGET_INLINE double      mcmath_rint<double>      (const double& x)      { return mc_rint  (x);    }
template <>        MC_TARGET_INLINE long double mcmath_rint<long double> (const long double& x) { return mc_rintl (x);    }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_rint (const float x)       { return mc_rintf (x); }
MC_TARGET_ALIAS double      mcmath_rint (const double x)      { return mc_rint  (x); }
MC_TARGET_ALIAS long double mcmath_rint (const long double x) { return mc_rintl (x); }
#	elif MC_TARGET_C11
#	define mcmath_rint(x) _Generic(x \
	, float       : mc_rintf         \
	, double      : mc_rint          \
	, long double : mc_rintl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_rint(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_rintf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_rint  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_rintl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_rint(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_rintf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_rint  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_rintl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_lrint -

#	ifndef mcmath_lrint
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE long mcmath_lrint              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE long mcmath_lrint<float>       (const float& x)       { return mc_lrintf (x);   }
template <>        MC_TARGET_INLINE long mcmath_lrint<double>      (const double& x)      { return mc_lrint  (x);   }
template <>        MC_TARGET_INLINE long mcmath_lrint<long double> (const long double& x) { return mc_lrintl (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS long mcmath_lrint (const float x)       { return mc_lrintf (x); }
MC_TARGET_ALIAS long mcmath_lrint (const double x)      { return mc_lrint  (x); }
MC_TARGET_ALIAS long mcmath_lrint (const long double x) { return mc_lrintl (x); }
#	elif MC_TARGET_C11
#	define mcmath_lrint(x) _Generic(x \
	, float       : mc_lrintf         \
	, double      : mc_lrint          \
	, long double : mc_lrintl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_lrint(x) \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_lrintf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_lrint  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_lrintl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_lrint(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_lrintf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_lrint  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_lrintl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_llrint -

#	ifndef mcmath_llrint
#	if MC_TARGET_CPP11
template <class T> MC_TARGET_INLINE long long mcmath_llrint              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE long long mcmath_llrint<float>       (const float& x)       { return mc_llrintf (x);  }
template <>        MC_TARGET_INLINE long long mcmath_llrint<double>      (const double& x)      { return mc_llrint  (x);  }
template <>        MC_TARGET_INLINE long long mcmath_llrint<long double> (const long double& x) { return mc_llrintl (x);  }
#	elif MC_TARGET_C99 && MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS long long mcmath_llrint (const float x)       { return mc_llrintf (x); }
MC_TARGET_ALIAS long long mcmath_llrint (const double x)      { return mc_llrint  (x); }
MC_TARGET_ALIAS long long mcmath_llrint (const long double x) { return mc_llrintl (x); }
#	elif MC_TARGET_C11
#	define mcmath_llrint(x) _Generic(x \
	, float       : mc_llrintf         \
	, double      : mc_llrint          \
	, long double : mc_llrintl         \
) (x)
#	elif MC_TARGET_C99 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_llrint(x) \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_llrintf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_llrint  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_llrintl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	elif MC_TARGET_C99
#	define mcmath_llrint(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_llrintf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_llrint  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_llrintl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	else
#	define mcmath_llrint mcmath_lrint
#	endif
#	endif

#pragma mark - mcmath_round -

#	ifndef mcmath_round
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_round              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_round<float>       (const float& x)       { return mc_roundf (x);   }
template <>        MC_TARGET_INLINE double      mcmath_round<double>      (const double& x)      { return mc_round  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_round<long double> (const long double& x) { return mc_roundl (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_round (const float x)       { return mc_roundf (x); }
MC_TARGET_ALIAS double      mcmath_round (const double x)      { return mc_round  (x); }
MC_TARGET_ALIAS long double mcmath_round (const long double x) { return mc_roundl (x); }
#	elif MC_TARGET_C11
#	define mcmath_round(x) _Generic(x \
	, float       : mc_roundf         \
	, double      : mc_round          \
	, long double : mc_roundl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_round(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_roundf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_round  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_roundl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_round(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_roundf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_round  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_roundl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_lround -

#	ifndef mcmath_lround
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE long mcmath_lround              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE long mcmath_lround<float>       (const float& x)       { return mc_lroundf (x);  }
template <>        MC_TARGET_INLINE long mcmath_lround<double>      (const double& x)      { return mc_lround  (x);  }
template <>        MC_TARGET_INLINE long mcmath_lround<long double> (const long double& x) { return mc_lroundl (x);  }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS long mcmath_lround (const float x)       { return mc_lroundf (x); }
MC_TARGET_ALIAS long mcmath_lround (const double x)      { return mc_lround  (x); }
MC_TARGET_ALIAS long mcmath_lround (const long double x) { return mc_lroundl (x); }
#	elif MC_TARGET_C11
#	define mcmath_lround(x) _Generic(x \
	, float       : mc_lroundf         \
	, double      : mc_lround          \
	, long double : mc_lroundl         \
) (x)
#	else
#	define mcmath_lround(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_lroundf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_lround  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_lroundl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_llround -

#	ifndef mcmath_llround
#	if MC_TARGET_CPP11
template <class T> MC_TARGET_INLINE long long mcmath_llround              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE long long mcmath_llround<float>       (const float& x)       { return mc_llroundf (x); }
template <>        MC_TARGET_INLINE long long mcmath_llround<double>      (const double& x)      { return mc_llround  (x); }
template <>        MC_TARGET_INLINE long long mcmath_llround<long double> (const long double& x) { return mc_llroundl (x); }
#	elif MC_TARGET_C99 && MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS long long mcmath_llround (const float x)       { return mc_llroundf (x); }
MC_TARGET_ALIAS long long mcmath_llround (const double x)      { return mc_llround  (x); }
MC_TARGET_ALIAS long long mcmath_llround (const long double x) { return mc_llroundl (x); }
#	elif MC_TARGET_C11
#	define mcmath_llround(x) _Generic(x \
	, float       : mc_lroundf          \
	, double      : mc_lround           \
	, long double : mc_lroundl          \
) (x)
#	elif MC_TARGET_C99 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_llrint(x) \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_llroundf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_llround  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_llroundl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	elif MC_TARGET_C99
#	define mcmath_llround(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_llroundf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_llround  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_llroundl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	else
#	define mcmath_llround(x) mcmath_lround(x)
#	endif
#	endif

#pragma mark - mcmath_trunc -

#	ifndef mcmath_trunc
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_trunc              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_trunc<float>       (const float& x)       { return mc_truncf (x);   }
template <>        MC_TARGET_INLINE double      mcmath_trunc<double>      (const double& x)      { return mc_trunc  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_trunc<long double> (const long double& x) { return mc_truncl (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_trunc (const float x)       { return mc_truncf (x); }
MC_TARGET_ALIAS double      mcmath_trunc (const double x)      { return mc_trunc  (x); }
MC_TARGET_ALIAS long double mcmath_trunc (const long double x) { return mc_truncl (x); }
#	elif MC_TARGET_C11
#	define mcmath_trunc(x) _Generic(x \
	, float       : mc_truncf         \
	, double      : mc_trunc          \
	, long double : mc_truncl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_trunc(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_truncf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_trunc  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_truncl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_trunc(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_truncf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_trunc  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_truncl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#	ifndef mcmath_fmod
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_fmod              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_fmod<float>       (const float& x, const float& y)             { return mc_fmodf (x, y);               }
template <>        MC_TARGET_INLINE double      mcmath_fmod<double>      (const double& x, const double& y)           { return mc_fmod  (x, y);               }
template <>        MC_TARGET_INLINE long double mcmath_fmod<long double> (const long double& x, const long double& y) { return mc_fmodl (x, y);               }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_fmod (const float x, const float y)             { return mc_fmodf (x, y); }
MC_TARGET_ALIAS double      mcmath_fmod (const double x, const double y)           { return mc_fmod  (x, y); }
MC_TARGET_ALIAS long double mcmath_fmod (const long double x, const long double y) { return mc_fmodl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_fmod(x, y) _Generic(x \
	, float       : mc_fmodf            \
	, double      : mc_fmod             \
	, long double : mc_fmodl            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_fmod(x, y) _Generic((x)+(y) \
	, float       : mc_fmodf                  \
	, double      : mc_fmod                   \
	, long double : mc_fmodl                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_fmod(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_fmodf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_fmod  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_fmodl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_fmod(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_fmodf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_fmod  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_fmodl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#	ifndef mcmath_remainder
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_remainder              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_remainder<float>       (const float& x, const float& y)             { return mc_remainderf (x, y);               }
template <>        MC_TARGET_INLINE double      mcmath_remainder<double>      (const double& x, const double& y)           { return mc_remainder  (x, y);               }
template <>        MC_TARGET_INLINE long double mcmath_remainder<long double> (const long double& x, const long double& y) { return mc_remainderl (x, y);               }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_remainder (const float x, const float y)             { return mc_remainderf (x, y); }
MC_TARGET_ALIAS double      mcmath_remainder (const double x, const double y)           { return mc_remainder  (x, y); }
MC_TARGET_ALIAS long double mcmath_remainder (const long double x, const long double y) { return mc_remainderl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_remainder(x, y) _Generic(x \
	, float       : mc_remainderf            \
	, double      : mc_remainder             \
	, long double : mc_remainderl            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_remainder(x, y) _Generic((x)+(y) \
	, float       : mc_remainderf                  \
	, double      : mc_remainder                   \
	, long double : mc_remainderl                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_remainder(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_remainderf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_remainder  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_remainderl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_remainder(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_remainderf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_remainder  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_remainderl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_remquo -

#	ifndef mcmath_remquo
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_remquo              (const T& x, const T& y, int * q)                     { mc_unused(x); mc_unused(y); mc_unused(q); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_remquo<float>       (const float& x, const float& y, int * q)             { return mc_remquof (x, y, q);                        }
template <>        MC_TARGET_INLINE double      mcmath_remquo<double>      (const double& x, const double& y, int * q)           { return mc_remquo  (x, y, q);                        }
template <>        MC_TARGET_INLINE long double mcmath_remquo<long double> (const long double& x, const long double& y, int * q) { return mc_remquol (x, y, q);                        }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_remquo (const float x, const float y, int * q)             { return mc_remquof (x, y, q); }
MC_TARGET_ALIAS double      mcmath_remquo (const double x, const double y, int * q)           { return mc_remquo  (x, y, q); }
MC_TARGET_ALIAS long double mcmath_remquo (const long double x, const long double y, int * q) { return mc_remquol (x, y, q); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_remquo(x, y, q) _Generic(x \
	, float       : mc_remquof               \
	, double      : mc_remquo                \
	, long double : mc_remquol               \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y), mc_cast_expr(int *, q))
#	elif MC_TARGET_C11
#	define mcmath_remquo(x, y, q) _Generic((x)+(y) \
	, float       : mc_remquof                     \
	, double      : mc_remquo                      \
	, long double : mc_remquol                     \
) ((x), (y), mc_cast_expr(int *, q))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_remquo(x, y, q) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_remquof (mc_cast_expr(const float, x), mc_cast_expr(const float, y), mc_cast_expr(int *, q))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_remquo  (mc_cast_expr(const double, x), mc_cast_expr(const double, y), mc_cast_expr(int *, q))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_remquol (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y), mc_cast_expr(int *, q)) \
		: 0 \
	))
#	else
#	define mcmath_remquo(x, y, q) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_remquof (mc_cast_expr(const float, x), mc_cast_expr(const float, y), mc_cast_expr(int *, q))             \
		: sizeof(x) == sizeof(double)      ? mc_remquo  (mc_cast_expr(const double, x), mc_cast_expr(const double, y), mc_cast_expr(int *, q))           \
		: sizeof(x) == sizeof(long double) ? mc_remquol (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y), mc_cast_expr(int *, q)) \
		: 0 \
	)
#	endif
#	endif

#	ifndef mcmath_copysign
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_copysign              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_copysign<float>       (const float& x, const float& y)             { return mc_copysignf (x, y);           }
template <>        MC_TARGET_INLINE double      mcmath_copysign<double>      (const double& x, const double& y)           { return mc_copysign  (x, y);           }
template <>        MC_TARGET_INLINE long double mcmath_copysign<long double> (const long double& x, const long double& y) { return mc_copysignl (x, y);           }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_copysign (const float x, const float y)             { return mc_copysignf (x, y); }
MC_TARGET_ALIAS double      mcmath_copysign (const double x, const double y)           { return mc_copysign  (x, y); }
MC_TARGET_ALIAS long double mcmath_copysign (const long double x, const long double y) { return mc_copysignl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_copysign(x, y) _Generic(x \
	, float       : mc_copysignf            \
	, double      : mc_copysign             \
	, long double : mc_copysignl            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_copysign(x, y) _Generic((x)+(y) \
	, float       : mc_copysignf                  \
	, double      : mc_copysign                   \
	, long double : mc_copysignl                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_copysign(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_copysignf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_copysign  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_copysignl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_copysign(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_copysignf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_copysign  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_copysignl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_nan -

#	ifndef mcmath_nan
#	define mcmath_nan(x, p) \
	(x) = ( \
		  sizeof(x) == sizeof(float)       ? nanf (mc_cast_expr(const char *, p)) \
		: sizeof(x) == sizeof(double)      ? nan  (mc_cast_expr(const char *, p)) \
		: sizeof(x) == sizeof(long double) ? nanl (mc_cast_expr(const char *, p)) \
		: 0 \
	)
#	endif

#	ifndef mcmath_nextafter
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_nextafter              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_nextafter<float>       (const float& x, const float& y)             { return mc_nextafterf (x, y);          }
template <>        MC_TARGET_INLINE double      mcmath_nextafter<double>      (const double& x, const double& y)           { return mc_nextafter  (x, y);          }
template <>        MC_TARGET_INLINE long double mcmath_nextafter<long double> (const long double& x, const long double& y) { return mc_nextafterl (x, y);          }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_nextafter (const float x, const float y)             { return mc_nextafterf (x, y); }
MC_TARGET_ALIAS double      mcmath_nextafter (const double x, const double y)           { return mc_nextafter  (x, y); }
MC_TARGET_ALIAS long double mcmath_nextafter (const long double x, const long double y) { return mc_nextafterl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_nextafter(x, y) _Generic(x \
	, float       : mc_nextafterf            \
	, double      : mc_nextafter             \
	, long double : mc_nextafterl            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_nextafter(x, y) _Generic((x)+(y) \
	, float       : mc_nextafterf                  \
	, double      : mc_nextafter                   \
	, long double : mc_nextafterl                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_nextafter(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_nextafterf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_nextafter  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_nextafterl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_nextafter(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_nextafterf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_nextafter  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_nextafterl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#	ifndef mcmath_nexttoward
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_nexttoward              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_nexttoward<float>       (const float& x, const float& y)             { return mc_nexttowardf (x, y);         }
template <>        MC_TARGET_INLINE double      mcmath_nexttoward<double>      (const double& x, const double& y)           { return mc_nexttoward  (x, y);         }
template <>        MC_TARGET_INLINE long double mcmath_nexttoward<long double> (const long double& x, const long double& y) { return mc_nexttowardl (x, y);         }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_nexttoward (const float x, const float y)             { return mc_nexttowardf (x, y); }
MC_TARGET_ALIAS double      mcmath_nexttoward (const double x, const double y)           { return mc_nexttoward  (x, y); }
MC_TARGET_ALIAS long double mcmath_nexttoward (const long double x, const long double y) { return mc_nexttowardl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_nexttoward(x, y) _Generic(x \
	, float       : mc_nexttowardf            \
	, double      : mc_nexttoward             \
	, long double : mc_nexttowardl            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_nexttoward(x, y) _Generic((x)+(y) \
	, float       : mc_nexttowardf                  \
	, double      : mc_nexttoward                   \
	, long double : mc_nexttowardl                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_nexttoward(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_nexttowardf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_nexttoward  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_nexttowardl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_nexttoward(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_nexttowardf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_nexttoward  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_nexttowardl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#	ifndef mcmath_fdim
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_fdim              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_fdim<float>       (const float& x, const float& y)             { return mc_fdimf (x, y);               }
template <>        MC_TARGET_INLINE double      mcmath_fdim<double>      (const double& x, const double& y)           { return mc_fdim  (x, y);               }
template <>        MC_TARGET_INLINE long double mcmath_fdim<long double> (const long double& x, const long double& y) { return mc_fdiml (x, y);               }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_fdim (const float x, const float y)             { return mc_fdimf (x, y); }
MC_TARGET_ALIAS double      mcmath_fdim (const double x, const double y)           { return mc_fdim  (x, y); }
MC_TARGET_ALIAS long double mcmath_fdim (const long double x, const long double y) { return mc_fdiml (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_fdim(x, y) _Generic(x \
	, float       : mc_fdimf            \
	, double      : mc_fdim             \
	, long double : mc_fdiml            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_fdim(x, y) _Generic((x)+(y) \
	, float       : mc_fdimf                  \
	, double      : mc_fdim                   \
	, long double : mc_fdiml                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_fdim(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_fdimf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_fdim  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_fdiml (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_fdim(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_fdimf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_fdim  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_fdiml (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_fma -

#	ifndef mcmath_fma
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_fma              (const T& x, const T& y, const T& z)                               { mc_unused(x); mc_unused(y); mc_unused(z); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_fma<float>       (const float& x, const float& y, const float& z)                   { return mc_fmaf (x, y, z);                           }
template <>        MC_TARGET_INLINE double      mcmath_fma<double>      (const double& x, const double& y, const double& z)                { return mc_fma  (x, y, z);                           }
template <>        MC_TARGET_INLINE long double mcmath_fma<long double> (const long double& x, const long double& y, const long double& z) { return mc_fmal (x, y, z);                           }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_fma (const float x, const float y, const float z)                   { return mc_fmaf (x, y, z); }
MC_TARGET_ALIAS double      mcmath_fma (const double x, const double y, const double z)                { return mc_fma  (x, y, z); }
MC_TARGET_ALIAS long double mcmath_fma (const long double x, const long double y, const long double z) { return mc_fmal (x, y, z); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_fma(x, y, z) _Generic(x \
	, float       : mc_fmaf               \
	, double      : mc_fma                \
	, long double : mc_fmal               \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y), mc_cast_expr(MC_TARGET_TYPEOF(x), z))
#	elif MC_TARGET_C11
#	define mcmath_fma(x, y, z) _Generic((x)+(y)+(z) \
	, float       : mc_fmaf                         \
	, double      : mc_fma                          \
	, long double : mc_fmal                         \
) ((x), (y), (z))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_fma(x, y, z) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_fmaf (mc_cast_expr(const float, x), mc_cast_expr(const float, y), mc_cast_expr(const float, z))                   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_fma  (mc_cast_expr(const double, x), mc_cast_expr(const double, y), mc_cast_expr(const double, z))                \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_fmal (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y), mc_cast_expr(const long double, z)) \
		: 0 \
	))
#	else
#	define mcmath_fma(x, y, z) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_fmaf (mc_cast_expr(const float, x), mc_cast_expr(const float, y), mc_cast_expr(const float, z))                   \
		: sizeof(x) == sizeof(double)      ? mc_fma  (mc_cast_expr(const double, x), mc_cast_expr(const double, y), mc_cast_expr(const double, z))                \
		: sizeof(x) == sizeof(long double) ? mc_fmal (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y), mc_cast_expr(const long double, z)) \
		: 0 \
	)
#	endif
#	endif

#	ifndef mcmath_fmax
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_fmax              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_fmax<float>       (const float& x, const float& y)             { return mc_fmaxf (x, y);               }
template <>        MC_TARGET_INLINE double      mcmath_fmax<double>      (const double& x, const double& y)           { return mc_fmax  (x, y);               }
template <>        MC_TARGET_INLINE long double mcmath_fmax<long double> (const long double& x, const long double& y) { return mc_fmaxl (x, y);               }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_fmax (const float x, const float y)             { return mc_fmaxf (x, y); }
MC_TARGET_ALIAS double      mcmath_fmax (const double x, const double y)           { return mc_fmax  (x, y); }
MC_TARGET_ALIAS long double mcmath_fmax (const long double x, const long double y) { return mc_fmaxl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_fmax(x, y) _Generic(x \
	, float       : mc_fmaxf            \
	, double      : mc_fmax             \
	, long double : mc_fmaxl            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_fmax(x, y) _Generic((x)+(y) \
	, float       : mc_fmaxf                  \
	, double      : mc_fmax                   \
	, long double : mc_fmaxl                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_fmax(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_fmaxf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_fmax  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_fmaxl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_fmax(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_fmaxf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_fmax  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_fmaxl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_fmin -

#	ifndef mcmath_fmin
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_fmin              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_fmin<float>       (const float& x, const float& y)             { return mc_fminf (x, y);               }
template <>        MC_TARGET_INLINE double      mcmath_fmin<double>      (const double& x, const double& y)           { return mc_fmin  (x, y);               }
template <>        MC_TARGET_INLINE long double mcmath_fmin<long double> (const long double& x, const long double& y) { return mc_fminl (x, y);               }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_fmin (const float x, const float y)             { return mc_fminf (x, y); }
MC_TARGET_ALIAS double      mcmath_fmin (const double x, const double y)           { return mc_fmin  (x, y); }
MC_TARGET_ALIAS long double mcmath_fmin (const long double x, const long double y) { return mc_fminl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_fmin(x, y) _Generic(x \
	, float       : mc_fminf            \
	, double      : mc_fmin             \
	, long double : mc_fminl            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_fmin(x, y) _Generic((x)+(y) \
	, float       : mc_fminf                  \
	, double      : mc_fmin                   \
	, long double : mc_fminl                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_fmin(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_fminf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_fmin  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_fminl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_fmin(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_fminf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_fmin  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_fminl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_abs -

#	ifndef mcmath_abs
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T                  mcmath_abs                      (const T& x)                  { return ::std::abs(x); }
template <>        MC_TARGET_INLINE float              mcmath_abs<float>               (const float& x)              { return mc_fabsf  (x); }
template <>        MC_TARGET_INLINE double             mcmath_abs<double>              (const double& x)             { return mc_fabs   (x); }
template <>        MC_TARGET_INLINE long double        mcmath_abs<long double>         (const long double& x)        { return mc_fabsl  (x); }

template <>        MC_TARGET_INLINE signed char        mcmath_abs<signed char>         (const signed char& x)        { return mc_babs   (x); }
template <>        MC_TARGET_INLINE short              mcmath_abs<short>               (const short& x)              { return mc_sabs   (x); }
template <>        MC_TARGET_INLINE int                mcmath_abs<int>                 (const int& x)                { return mc_iabs   (x); }
template <>        MC_TARGET_INLINE long               mcmath_abs<long>                (const long& x)               { return mc_labs   (x); }

template <>        MC_TARGET_INLINE unsigned char      mcmath_abs<unsigned char>       (const unsigned char& x)      { return x;             }
template <>        MC_TARGET_INLINE unsigned short     mcmath_abs<unsigned short>      (const unsigned short& x)     { return x;             }
template <>        MC_TARGET_INLINE unsigned int       mcmath_abs<unsigned int>        (const unsigned int& x)       { return x;             }
template <>        MC_TARGET_INLINE unsigned long      mcmath_abs<unsigned long>       (const unsigned long& x)      { return x;             }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE long long          mcmath_abs<long long>           (const long long& x)          { return x;             }
template <>        MC_TARGET_INLINE unsigned long long mcmath_abs<unsigned  long long> (const unsigned long long& x) { return x;             }
#	endif
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float              mcmath_abs (const float x)              { return mc_fabsf (x); }
MC_TARGET_ALIAS double             mcmath_abs (const double x)             { return mc_fabs  (x); }
MC_TARGET_ALIAS long double        mcmath_abs (const long double x)        { return mc_fabsl (x); }
MC_TARGET_ALIAS signed char        mcmath_abs (const signed char x)        { return mc_babs  (x); }
MC_TARGET_ALIAS short              mcmath_abs (const short x)              { return mc_sabs  (x); }
MC_TARGET_ALIAS int                mcmath_abs (const int x)                { return mc_iabs  (x); }
MC_TARGET_ALIAS long               mcmath_abs (const long x)               { return mc_labs  (x); }
MC_TARGET_ALIAS unsigned char      mcmath_abs (const unsigned char x)      { return x;            }
MC_TARGET_ALIAS unsigned short     mcmath_abs (const unsigned short x)     { return x;            }
MC_TARGET_ALIAS unsigned int       mcmath_abs (const unsigned int x)       { return x;            }
MC_TARGET_ALIAS unsigned long      mcmath_abs (const unsigned long x)      { return x;            }
#	if MC_TARGET_C99
MC_TARGET_ALIAS long long          mcmath_abs (const long long x)          { return mc_llabs (x); }
MC_TARGET_ALIAS unsigned long long mcmath_abs (const unsigned long long x) { return x;            }
#	endif
#	elif MC_TARGET_C11
#	define mcmath_abs(x) _Generic(x            \
	, float              : mc_fabsf            \
	, double             : mc_fabs             \
	, long double        : mc_fabsl            \
	, signed char        : mc_babs             \
	, short              : mc_sabs             \
	, int                : mc_iabs             \
	, long               : mc_labs             \
	, long long          : mc_llabs            \
	, unsigned char      : MC_TARGET_UCHAR     \
	, unsigned short     : MC_TARGET_USHORT    \
	, unsigned int       : MC_TARGET_UINT      \
	, unsigned long      : MC_TARGET_ULONG     \
	, unsigned long long : MC_TARGET_ULONGLONG \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#		if MC_TARGET_C99
#		define mcmath_abs(x) mc_cast(MC_TARGET_TYPEOF(x), \
		( \
			  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? mc_fabsf (mc_cast_expr(const float, x))       \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? mc_fabs  (mc_cast_expr(const double, x))      \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? mc_fabsl (mc_cast_expr(const long double, x)) \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? mc_babs  (mc_cast_expr(signed char, x))       \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? mc_sabs  (mc_cast_expr(short, x))             \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? mc_iabs  (mc_cast_expr(int,x))                \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? mc_labs  (mc_cast_expr(long, x))              \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long long)          ? mc_llabs (mc_cast_expr(long long, x))         \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? (x)                                           \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? (x)                                           \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? (x)                                           \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? (x)                                           \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long long) ? (x)                                           \
			: 0 \
		))
#		else
#		define mcmath_abs(x) mc_cast(MC_TARGET_TYPEOF(x), \
		( \
			  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? mc_fabsf (mc_cast_expr(const float, x))       \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? mc_fabs  (mc_cast_expr(const double, x))      \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? mc_fabsl (mc_cast_expr(const long double, x)) \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? mc_babs  (mc_cast_expr(signed char, x))       \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? mc_sabs  (mc_cast_expr(short, x))             \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? mc_iabs  (mc_cast_expr(int,x))                \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? mc_labs  (mc_cast_expr(long, x))              \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? (x)                                           \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? (x)                                           \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? (x)                                           \
			: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? (x)                                           \
			: 0 \
		))
#		endif
#	else
#	define mcmath_abs(x) mc_absmag(x)
#	endif
#	endif

#pragma mark - mcmath_max -

#	ifndef mcmath_max
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_max              (const T& x, const T& y)                     { return ::std::max (x, y); }
template <>        MC_TARGET_INLINE float       mcmath_max<float>       (const float& x, const float& y)             { return mc_fmaxf   (x, y); }
template <>        MC_TARGET_INLINE double      mcmath_max<double>      (const double& x, const double& y)           { return mc_fmax    (x, y); }
template <>        MC_TARGET_INLINE long double mcmath_max<long double> (const long double& x, const long double& y) { return mc_fmaxl   (x, y); }
#	elif MC_TARGET_C99 && MC_TARGET_HAVE_AUTOTYPE
#	define mcmath_max(a, b) \
	__extension__ ({ \
		MC_TARGET_AUTOTYPE __mcmath_max_aa = (a);  \
		MC_TARGET_AUTOTYPE __mcmath_max_bb = (b);  \
		mc_cast(MC_TARGET_TYPEOF(__mcmath_max_aa), \
		( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF((__mcmath_max_aa)), float)       ? mc_fmaxf (mc_cast(float, __mcmath_max_aa)       , mc_cast(float, __mcmath_max_bb))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF((__mcmath_max_aa)), double)      ? mc_fmax  (mc_cast(double, __mcmath_max_aa)      , mc_cast(double, __mcmath_max_bb))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF((__mcmath_max_aa)), long double) ? mc_fmaxl (mc_cast(long double, __mcmath_max_aa) , mc_cast(long double, __mcmath_max_bb)) \
		: (__mcmath_max_aa > __mcmath_max_bb ? __mcmath_max_aa : __mcmath_max_bb)                                                                                         \
		)); \
	})
#	elif MC_TARGET_C99 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_max(a, b) \
	__extension__ ({ \
		MC_TARGET_TYPEOF((a)) __mcmath_max_aa = (a); \
		MC_TARGET_TYPEOF((b)) __mcmath_max_bb = (b); \
		mc_cast(MC_TARGET_TYPEOF(__mcmath_max_aa),   \
		( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF((__mcmath_max_aa)), float)       ? mc_fmaxf (mc_cast(float, __mcmath_max_aa)       , mc_cast(float, __mcmath_max_bb))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF((__mcmath_max_aa)), double)      ? mc_fmax  (mc_cast(double, __mcmath_max_aa)      , mc_cast(double, __mcmath_max_bb))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF((__mcmath_max_aa)), long double) ? mc_fmaxl (mc_cast(long double, __mcmath_max_aa) , mc_cast(long double, __mcmath_max_bb)) \
		: (__mcmath_max_aa > __mcmath_max_bb ? __mcmath_max_aa : __mcmath_max_bb)                                                                                         \
		)); \
	})
#	else
#	define mcmath_max(a, b) mc_maxmag(a, b)
#	endif
#	endif

#pragma mark - mcmath_min -

#	ifndef mcmath_min
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_min              (const T& x, const T& y)                     { return ::std::min (x, y); }
template <>        MC_TARGET_INLINE float       mcmath_min<float>       (const float& x, const float& y)             { return mc_fminf   (x, y); }
template <>        MC_TARGET_INLINE double      mcmath_min<double>      (const double& x, const double& y)           { return mc_fmin    (x, y); }
template <>        MC_TARGET_INLINE long double mcmath_min<long double> (const long double& x, const long double& y) { return mc_fminl   (x, y); }
#	elif MC_TARGET_C99 && MC_TARGET_HAVE_AUTOTYPE
#	define mcmath_min(a, b) \
	__extension__ ({ \
		MC_TARGET_AUTOTYPE __mcmath_min_aa = (a);  \
		MC_TARGET_AUTOTYPE __mcmath_min_bb = (b);  \
		mc_cast(MC_TARGET_TYPEOF(__mcmath_min_aa), \
		( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF((__mcmath_min_aa)), float)       ? mc_fminf (mc_cast(float, __mcmath_min_aa)       , mc_cast(float, __mcmath_min_bb))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF((__mcmath_min_aa)), double)      ? mc_fmin  (mc_cast(double, __mcmath_min_aa)      , mc_cast(double, __mcmath_min_bb))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF((__mcmath_min_aa)), long double) ? mc_fminl (mc_cast(long double, __mcmath_min_aa) , mc_cast(long double, __mcmath_min_bb)) \
		: (__mcmath_min_aa < __mcmath_min_bb ? __mcmath_min_aa : __mcmath_min_bb)                                                                                         \
		)); \
	})
#	elif MC_TARGET_C99 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_min(a, b) \
	__extension__ ({ \
		MC_TARGET_TYPEOF((a)) __mcmath_min_aa = (a); \
		MC_TARGET_TYPEOF((b)) __mcmath_min_bb = (b); \
		mc_cast(MC_TARGET_TYPEOF(__mcmath_min_aa), \
		( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF((__mcmath_min_aa)), float)       ? mc_fminf (mc_cast(float, __mcmath_min_aa), mc_cast(float, __mcmath_min_bb))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF((__mcmath_min_aa)), double)      ? mc_fmin  (mc_cast(double, __mcmath_min_aa), mc_cast(double, __mcmath_min_bb))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF((__mcmath_min_aa)), long double) ? mc_fminl (mc_cast(long double, __mcmath_min_aa), mc_cast(long double, __mcmath_min_bb)) \
		: (__mcmath_min_aa < __mcmath_min_bb ? __mcmath_min_aa : __mcmath_min_bb)                                                                                        \
		)); \
	})
#	else
#	define mcmath_min(a, b) mc_minmag(a, b)
#	endif
#	endif

#	if MC_TARGET_CPP11
#		define mcmath_isgreater(x, y)      ::std::isgreater(x, y)
#		define mcmath_isgreaterequal(x, y) ::std::isgreaterequal(x, y)
#		define mcmath_isless(x, y)         ::std::isless(x, y)
#		define mcmath_islessequal(x, y)    ::std::islessequal(x, y)
#		define mcmath_islessgreater(x, y)  ::std::islessgreater(x, y)
#		define mcmath_isunordered(x, y)    ::std::isunordered(x, y)
#	elif MC_TARGET_CPP98
#		define mcmath_isgreater(x, y)      ::isgreater(x, y)
#		define mcmath_isgreaterequal(x, y) ::isgreaterequal(x, y)
#		define mcmath_isless(x, y)         ::isless(x, y)
#		define mcmath_islessequal(x, y)    ::islessequal(x, y)
#		define mcmath_islessgreater(x, y)  ::islessgreater(x, y)
#		define mcmath_isunordered(x, y)    ::isunordered(x, y)
#	else
#		define mcmath_isgreater(x, y)      isgreater(x, y)
#		define mcmath_isgreaterequal(x, y) isgreaterequal(x, y)
#		define mcmath_isless(x, y)         isless(x, y)
#		define mcmath_islessequal(x, y)    islessequal(x, y)
#		define mcmath_islessgreater(x, y)  islessgreater(x, y)
#		define mcmath_isunordered(x, y)    isunordered(x, y)
#	endif

#	ifndef mcmath_i0
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_i0              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_i0<float>       (const float& x)       { return mc_i0f (x);      }
template <>        MC_TARGET_INLINE double      mcmath_i0<double>      (const double& x)      { return mc_i0  (x);      }
template <>        MC_TARGET_INLINE long double mcmath_i0<long double> (const long double& x) { return mc_i0l (x);      }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_i0 (const float x)       { return mc_i0f (x); }
MC_TARGET_ALIAS double      mcmath_i0 (const double x)      { return mc_i0  (x); }
MC_TARGET_ALIAS long double mcmath_i0 (const long double x) { return mc_i0l (x); }
#	elif MC_TARGET_C11
#	define mcmath_i0(x) _Generic(x \
	, float       : mc_i0f         \
	, double      : mc_i0          \
	, long double : mc_i0l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_i0(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_i0f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_i0  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_i0l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_i0(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_i0f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_i0  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_i0l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_i1 -

#	ifndef mcmath_i1
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_i1              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_i1<float>       (const float& x)       { return mc_i1f (x);      }
template <>        MC_TARGET_INLINE double      mcmath_i1<double>      (const double& x)      { return mc_i1  (x);      }
template <>        MC_TARGET_INLINE long double mcmath_i1<long double> (const long double& x) { return mc_i1l (x);      }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_i1 (const float x)       { return mc_i1f (x); }
MC_TARGET_ALIAS double      mcmath_i1 (const double x)      { return mc_i1  (x); }
MC_TARGET_ALIAS long double mcmath_i1 (const long double x) { return mc_i1l (x); }
#	elif MC_TARGET_C11
#	define mcmath_i1(x) _Generic(x \
	, float       : mc_i1f         \
	, double      : mc_i1          \
	, long double : mc_i1l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_i1(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_i1f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_i1  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_i1l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_i1(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_i1f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_i1  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_i1l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_in -

#	ifndef mcmath_in
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_in              (const int n, const T& x)           { mc_unused(n); mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_in<float>       (const int n, const float& x)       { return mc_inf (n, x);                 }
template <>        MC_TARGET_INLINE double      mcmath_in<double>      (const int n, const double& x)      { return mc_in  (n, x);                 }
template <>        MC_TARGET_INLINE long double mcmath_in<long double> (const int n, const long double& x) { return mc_inl (n, x);                 }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_in (const int n, float x)       { return mc_inf (n, x); }
MC_TARGET_ALIAS double      mcmath_in (const int n, double x)      { return mc_in  (n, x); }
MC_TARGET_ALIAS long double mcmath_in (const int n, long double x) { return mc_inl (n, x); }
#	elif MC_TARGET_C11
#	define mcmath_in(n, x) _Generic(x \
	, float       : mc_inf            \
	, double      : mc_in             \
	, long double : mc_inl            \
) (mc_cast_expr(const int, n), x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_in(n, x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_inf (mc_cast_expr(const int, n), mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_in  (mc_cast_expr(const int, n), mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_inl (mc_cast_expr(const int, n), mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_in(n, x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_inf (mc_cast_expr(const int, n), mc_cast_expr(const float, x)))       \
		: sizeof(x) == sizeof(double)      ? mc_in  (mc_cast_expr(const int, n), mc_cast_expr(const double, x))       \
		: sizeof(x) == sizeof(long double) ? mc_inl (mc_cast_expr(const int, n), mc_cast_expr(const long double, x))) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_j0 -

#	ifndef mcmath_j0
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_j0              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_j0<float>       (const float& x)       { return mc_j0f (x);      }
template <>        MC_TARGET_INLINE double      mcmath_j0<double>      (const double& x)      { return mc_j0  (x);      }
template <>        MC_TARGET_INLINE long double mcmath_j0<long double> (const long double& x) { return mc_j0l (x);      }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_j0 (const float x)       { return mc_j0f (x); }
MC_TARGET_ALIAS double      mcmath_j0 (const double x)      { return mc_j0  (x); }
MC_TARGET_ALIAS long double mcmath_j0 (const long double x) { return mc_j0l (x); }
#	elif MC_TARGET_C11
#	define mcmath_j0(x) _Generic(x \
	, float       : mc_j0f \
	, double      : mc_j0  \
	, long double : mc_j0l \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_j0(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_j0f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_j0  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_j0l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_j0(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_j0f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_j0  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_j0l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_j1 -

#	ifndef mcmath_j1
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_j1              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_j1<float>       (const float& x)       { return mc_j1f (x);      }
template <>        MC_TARGET_INLINE double      mcmath_j1<double>      (const double& x)      { return mc_j1  (x);      }
template <>        MC_TARGET_INLINE long double mcmath_j1<long double> (const long double& x) { return mc_j1l (x);      }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_j1 (const float x)       { return mc_j1f (x); }
MC_TARGET_ALIAS double      mcmath_j1 (const double x)      { return mc_j1  (x); }
MC_TARGET_ALIAS long double mcmath_j1 (const long double x) { return mc_j1l (x); }
#	elif MC_TARGET_C11
#	define mcmath_j1(x) _Generic(x \
	, float       : mc_j1f         \
	, double      : mc_j1          \
	, long double : mc_j1l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_j1(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_j1f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_j1  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_j1l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_j1(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_j1f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_j1  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_j1l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_jn -

#	ifndef mcmath_jn
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_jn              (const int n, const T& x)           { mc_unused(n); mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_jn<float>       (const int n, const float& x)       { return mc_jnf (n, x);                 }
template <>        MC_TARGET_INLINE double      mcmath_jn<double>      (const int n, const double& x)      { return mc_jn  (n, x);                 }
template <>        MC_TARGET_INLINE long double mcmath_jn<long double> (const int n, const long double& x) { return mc_jnl (n, x);                 }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_jn (const int n, float x)       { return mc_jnf (n, x); }
MC_TARGET_ALIAS double      mcmath_jn (const int n, double x)      { return mc_jn  (n, x); }
MC_TARGET_ALIAS long double mcmath_jn (const int n, long double x) { return mc_jnl (n, x); }
#	elif MC_TARGET_C11
#	define mcmath_jn(n, x) _Generic(x \
	, float       : mc_jnf            \
	, double      : mc_jn             \
	, long double : mc_jnl            \
) (mc_cast_expr(const int, n), x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_jn(n, x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_jnf (mc_cast_expr(const int, n), mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_jn  (mc_cast_expr(const int, n), mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_jnl (mc_cast_expr(const int, n), mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_jn(n, x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_jnf (mc_cast_expr(const int, n), mc_cast_expr(const float, x)))       \
		: sizeof(x) == sizeof(double)      ? mc_jn  (mc_cast_expr(const int, n), mc_cast_expr(const double, x))       \
		: sizeof(x) == sizeof(long double) ? mc_jnl (mc_cast_expr(const int, n), mc_cast_expr(const long double, x))) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_y0 -

#	ifndef mcmath_y0
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_y0              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_y0<float>       (const float& x)       { return mc_y0f (x);      }
template <>        MC_TARGET_INLINE double      mcmath_y0<double>      (const double& x)      { return mc_y0  (x);      }
template <>        MC_TARGET_INLINE long double mcmath_y0<long double> (const long double& x) { return mc_y0l (x);      }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_y0 (const float x)       { return mc_y0f (x); }
MC_TARGET_ALIAS double      mcmath_y0 (const double x)      { return mc_y0  (x); }
MC_TARGET_ALIAS long double mcmath_y0 (const long double x) { return mc_y0l (x); }
#	elif MC_TARGET_C11
#	define mcmath_y0(x) _Generic(x \
	, float       : mc_y0f         \
	, double      : mc_y0          \
	, long double : mc_y0l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_y0(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_y0f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_y0  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_y0l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_y0(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_y0f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_y0  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_y0l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_y1 -

#	ifndef mcmath_y1
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_y1              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_y1<float>       (const float& x)       { return mc_y1f (x);      }
template <>        MC_TARGET_INLINE double      mcmath_y1<double>      (const double& x)      { return mc_y1  (x);      }
template <>        MC_TARGET_INLINE long double mcmath_y1<long double> (const long double& x) { return mc_y1l (x);      }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_y1 (const float x)       { return mc_y1f (x); }
MC_TARGET_ALIAS double      mcmath_y1 (const double x)      { return mc_y1  (x); }
MC_TARGET_ALIAS long double mcmath_y1 (const long double x) { return mc_y1l (x); }
#	elif MC_TARGET_C11
#	define mcmath_y1(x) _Generic(x \
	, float       : mc_y1f         \
	, double      : mc_y1          \
	, long double : mc_y1l         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_y1(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_y1f (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_y1  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_y1l (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_y1(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_y1f (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_y1  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_y1l (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_yn -

#	ifndef mcmath_yn
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_yn              (const int n, const T& x)           { mc_unused(n); mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_yn<float>       (const int n, const float& x)       { return mc_ynf (n, x);                 }
template <>        MC_TARGET_INLINE double      mcmath_yn<double>      (const int n, const double& x)      { return mc_yn  (n, x);                 }
template <>        MC_TARGET_INLINE long double mcmath_yn<long double> (const int n, const long double& x) { return mc_ynl (n, x);                 }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_yn (const int n, float x)       { return mc_ynf (n, x); }
MC_TARGET_ALIAS double      mcmath_yn (const int n, double x)      { return mc_yn  (n, x); }
MC_TARGET_ALIAS long double mcmath_yn (const int n, long double x) { return mc_ynl (n, x); }
#	elif MC_TARGET_C11
#	define mcmath_yn(n, x) _Generic(x \
	, float       : mc_ynf            \
	, double      : mc_yn             \
	, long double : mc_ynl            \
) (mc_cast_expr(const int, n), x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_yn(n, x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_ynf (mc_cast_expr(const int, n), mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_yn  (mc_cast_expr(const int, n), mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_ynl (mc_cast_expr(const int, n), mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_yn(n, x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_ynf (mc_cast_expr(const int, n), mc_cast_expr(const float, x)))       \
		: sizeof(x) == sizeof(double)      ? mc_yn  (mc_cast_expr(const int, n), mc_cast_expr(const double, x))       \
		: sizeof(x) == sizeof(long double) ? mc_ynl (mc_cast_expr(const int, n), mc_cast_expr(const long double, x))) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_scalb -

#	ifndef mcmath_scalb
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_scalb              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_scalb<float>       (const float& x, const float& y)             { return mc_scalbf (x, y);              }
template <>        MC_TARGET_INLINE double      mcmath_scalb<double>      (const double& x, const double& y)           { return mc_scalb  (x, y);              }
template <>        MC_TARGET_INLINE long double mcmath_scalb<long double> (const long double& x, const long double& y) { return mc_scalbl (x, y);              }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_scalb (const float x, const float y)             { return mc_scalbf (x, y); }
MC_TARGET_ALIAS double      mcmath_scalb (const double x, const double y)           { return mc_scalb  (x, y); }
MC_TARGET_ALIAS long double mcmath_scalb (const long double x, const long double y) { return mc_scalbl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_scalb(x, y) _Generic(x \
	, float       : mc_scalbf            \
	, double      : mc_scalb             \
	, long double : mc_scalbl            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_scalb(x, y) _Generic((x)+(y) \
	, float       : mc_scalbf                  \
	, double      : mc_scalb                   \
	, long double : mc_scalbl                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_scalb(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_scalbf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_scalb  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_scalbl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_scalb(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_scalbf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_scalb  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_scalbl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_ibeta -

#	ifndef mcmath_ibeta
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_ibeta              (const T& a, const T& b, const T& x)                               { mc_unused(a); mc_unused(b); mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_ibeta<float>       (const float& a, const float& b, const float& x)                   { return mc_ibetaf (a, b, x);                         }
template <>        MC_TARGET_INLINE double      mcmath_ibeta<double>      (const double& a, const double& b, const double& x)                { return mc_ibeta  (a, b, x);                         }
template <>        MC_TARGET_INLINE long double mcmath_ibeta<long double> (const long double& a, const long double& b, const long double& x) { return mc_ibetal (a, b, x);                         }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_ibeta (const float a, const float b, const float x)                   { return mc_ibetaf (a, b, x); }
MC_TARGET_ALIAS double      mcmath_ibeta (const double a, const double b, const double x)                { return mc_ibeta  (a, b, x); }
MC_TARGET_ALIAS long double mcmath_ibeta (const long double a, const long double b, const long double x) { return mc_ibetal (a, b, x); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_ibeta(a, b, x) _Generic(a \
	, float       : mc_ibetaf               \
	, double      : mc_ibeta                \
	, long double : mc_ibetal               \
) (a, mc_cast_expr(MC_TARGET_TYPEOF(a), b), mc_cast_expr(MC_TARGET_TYPEOF(a), x))
#	elif MC_TARGET_C11
#	define mcmath_ibeta(a, b, x) _Generic((a)+(b)+(x) \
	, float       : mc_ibetaf                         \
	, double      : mc_ibeta                          \
	, long double : mc_ibetal                         \
) ((a), (b), (x))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_ibeta(a, b, x) mc_cast(MC_TARGET_TYPEOF(a), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), float)       ? mc_ibetaf (mc_cast_expr(const float, a), mc_cast_expr(const float, b), mc_cast_expr(const float, x))                   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), double)      ? mc_ibeta  (mc_cast_expr(const double, a), mc_cast_expr(const double, b), mc_cast_expr(const double, x))                \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), long double) ? mc_ibetal (mc_cast_expr(const long double, a), mc_cast_expr(const long double, b), mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_ibeta(a, b, x) \
	( \
		  sizeof(a) == sizeof(float)       ? mc_ibetaf (mc_cast_expr(const float, a), mc_cast_expr(const float, b), mc_cast_expr(const float, x))                   \
		: sizeof(a) == sizeof(double)      ? mc_ibeta  (mc_cast_expr(const double, a), mc_cast_expr(const double, b), mc_cast_expr(const double, x))                \
		: sizeof(a) == sizeof(long double) ? mc_ibetal (mc_cast_expr(const long double, a), mc_cast_expr(const long double, b), mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_lbeta -

#	ifndef mcmath_lbeta
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_lbeta              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_lbeta<float>       (const float& x, const float& y)             { return mc_lbetaf (x, y);              }
template <>        MC_TARGET_INLINE double      mcmath_lbeta<double>      (const double& x, const double& y)           { return mc_lbeta  (x, y);              }
template <>        MC_TARGET_INLINE long double mcmath_lbeta<long double> (const long double& x, const long double& y) { return mc_lbetal (x, y);              }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_lbeta (const float x, const float y)             { return mc_lbetaf (x, y); }
MC_TARGET_ALIAS double      mcmath_lbeta (const double x, const double y)           { return mc_lbeta  (x, y); }
MC_TARGET_ALIAS long double mcmath_lbeta (const long double x, const long double y) { return mc_lbetal (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_lbeta(x, y) _Generic(x \
	, float       : mc_lbetaf            \
	, double      : mc_lbeta             \
	, long double : mc_lbetal            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_lbeta(x, y) _Generic((x)+(y) \
	, float       : mc_lbetaf                  \
	, double      : mc_lbeta                   \
	, long double : mc_lbetal                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_lbeta(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_lbetaf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_lbeta  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_lbetal (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_lbeta(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_lbetaf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_lbeta  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_lbetal (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_beta -

#	ifndef mcmath_beta
#	if MC_TARGET_CPP17 && MC_TARGET_HAVE_BETAFN
#	define mcmath_betaf ::std::betaf
#	define mcmath_beta  ::std::beta
#	define mcmath_betal ::std::betal
#	else
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_beta              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_beta<float>       (const float& x, const float& y)             { return mc_betaf (x, y);               }
template <>        MC_TARGET_INLINE double      mcmath_beta<double>      (const double& x, const double& y)           { return mc_beta  (x, y);               }
template <>        MC_TARGET_INLINE long double mcmath_beta<long double> (const long double& x, const long double& y) { return mc_betal (x, y);               }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_beta (const float x, const float y)             { return mc_betaf (x, y); }
MC_TARGET_ALIAS double      mcmath_beta (const double x, const double y)           { return mc_beta  (x, y); }
MC_TARGET_ALIAS long double mcmath_beta (const long double x, const long double y) { return mc_betal (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_beta(x, y) _Generic(x \
	, float       : mc_betaf            \
	, double      : mc_beta             \
	, long double : mc_betal            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_beta(x, y) _Generic((x)+(y) \
	, float       : mc_betaf                  \
	, double      : mc_beta                  \
	, long double : mc_betal                 \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_beta(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_betaf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_beta  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_betal (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_beta(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_betaf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_beta  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_betal (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif
#	endif

#	define mcmath_betainc(a, b, x) mcmath_ibeta(a, b, x)
#	define mcmath_betaln(x, y)     mcmath_lbeta(x, y)

#pragma mark - mcmath_xlogy -

#	ifndef mcmath_xlogy
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_xlogy              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_xlogy<float>       (const float& x, const float& y)             { return mc_xlogyf (x, y);              }
template <>        MC_TARGET_INLINE double      mcmath_xlogy<double>      (const double& x, const double& y)           { return mc_xlogy  (x, y);              }
template <>        MC_TARGET_INLINE long double mcmath_xlogy<long double> (const long double& x, const long double& y) { return mc_xlogyl (x, y);              }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_xlogy (const float x, const float y)             { return mc_xlogyf (x, y); }
MC_TARGET_ALIAS double      mcmath_xlogy (const double x, const double y)           { return mc_xlogy  (x, y); }
MC_TARGET_ALIAS long double mcmath_xlogy (const long double x, const long double y) { return mc_xlogyl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_xlogy(x, y) _Generic(x \
	, float       : mc_xlogyf            \
	, double      : mc_xlogy             \
	, long double : mc_xlogyl            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_xlogy(x, y) _Generic((x)+(y) \
	, float       : mc_xlogyf                  \
	, double      : mc_xlogy                   \
	, long double : mc_xlogyl                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_xlogy(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_xlogyf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_xlogy  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_xlogyl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_xlogy(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_xlogyf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_xlogy  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_xlogyl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_xlog1py -

#	ifndef mcmath_xlog1py
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_xlog1py              (const T& x, const T& y)                     { mc_unused(x); mc_unused(y); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_xlog1py<float>       (const float& x, const float& y)             { return mc_xlog1pyf (x, y);            }
template <>        MC_TARGET_INLINE double      mcmath_xlog1py<double>      (const double& x, const double& y)           { return mc_xlog1py  (x, y);            }
template <>        MC_TARGET_INLINE long double mcmath_xlog1py<long double> (const long double& x, const long double& y) { return mc_xlog1pyl (x, y);            }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_xlog1py (const float x, const float y)             { return mc_xlog1pyf (x, y); }
MC_TARGET_ALIAS double      mcmath_xlog1py (const double x, const double y)           { return mc_xlog1py  (x, y); }
MC_TARGET_ALIAS long double mcmath_xlog1py (const long double x, const long double y) { return mc_xlog1pyl (x, y); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_xlog1py(x, y) _Generic(x \
	, float       : mc_xlog1pyf            \
	, double      : mc_xlog1py             \
	, long double : mc_xlog1pyl            \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y))
#	elif MC_TARGET_C11
#	define mcmath_xlog1py(x, y) _Generic((x)+(y) \
	, float       : mc_xlog1pyf                  \
	, double      : mc_xlog1py                   \
	, long double : mc_xlog1pyl                  \
) ((x), (y))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_xlog1py(x, y) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_xlog1pyf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_xlog1py  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_xlog1pyl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	))
#	else
#	define mcmath_xlog1py(x, y) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_xlog1pyf (mc_cast_expr(const float, x), mc_cast_expr(const float, y))             \
		: sizeof(x) == sizeof(double)      ? mc_xlog1py  (mc_cast_expr(const double, x), mc_cast_expr(const double, y))           \
		: sizeof(x) == sizeof(long double) ? mc_xlog1pyl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_xlog1px -

#	ifndef mcmath_xlog1px
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_xlog1px              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_xlog1px<float>       (const float& x)       { return mc_xlog1pxf (x); }
template <>        MC_TARGET_INLINE double      mcmath_xlog1px<double>      (const double& x)      { return mc_xlog1px  (x); }
template <>        MC_TARGET_INLINE long double mcmath_xlog1px<long double> (const long double& x) { return mc_xlog1pxl (x); }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_xlog1px (const float x)       { return mc_xlog1pxf (x); }
MC_TARGET_ALIAS double      mcmath_xlog1px (const double x)      { return mc_xlog1px  (x); }
MC_TARGET_ALIAS long double mcmath_xlog1px (const long double x) { return mc_xlog1pxl (x); }
#	elif MC_TARGET_C11
#	define mcmath_xlog1px(x) _Generic(x \
	, float       : mc_xlog1pxf         \
	, double      : mc_xlog1px          \
	, long double : mc_xlog1pxl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_xlog1px(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_xlog1pxf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_xlog1px  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_xlog1pxl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_xlog1px(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_xlog1pxf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_xlog1px  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_xlog1pxl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_xlogx -

#	ifndef mcmath_xlogx
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_xlogx              (const T& x)           { mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_xlogx<float>       (const float& x)       { return mc_xlogxf (x);   }
template <>        MC_TARGET_INLINE double      mcmath_xlogx<double>      (const double& x)      { return mc_xlogx  (x);   }
template <>        MC_TARGET_INLINE long double mcmath_xlogx<long double> (const long double& x) { return mc_xlogxl (x);   }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_xlogx (const float x)       { return mc_xlogxf (x); }
MC_TARGET_ALIAS double      mcmath_xlogx (const double x)      { return mc_xlogx  (x); }
MC_TARGET_ALIAS long double mcmath_xlogx (const long double x) { return mc_xlogxl (x); }
#	elif MC_TARGET_C11
#	define mcmath_xlogx(x) _Generic(x \
	, float       : mc_xlogxf         \
	, double      : mc_xlogx          \
	, long double : mc_xlogxl         \
) (x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_xlogx(x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_xlogxf (mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_xlogx  (mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_xlogxl (mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_xlogx(x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_xlogxf (mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_xlogx  (mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_xlogxl (mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#	ifndef mcmath_logradix
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_logradix              (const unsigned int n, const T& x)           { mc_unused(n); mc_unused(x); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_logradix<float>       (const unsigned int n, const float& x)       { return mc_logradixf (n, x);           }
template <>        MC_TARGET_INLINE double      mcmath_logradix<double>      (const unsigned int n, const double& x)      { return mc_logradix  (n, x);           }
template <>        MC_TARGET_INLINE long double mcmath_logradix<long double> (const unsigned int n, const long double& x) { return mc_logradixl (n, x);           }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_logradix (const unsigned int n, float x)       { return mc_logradixf (n, x); }
MC_TARGET_ALIAS double      mcmath_logradix (const unsigned int n, double x)      { return mc_logradix  (n, x); }
MC_TARGET_ALIAS long double mcmath_logradix (const unsigned int n, long double x) { return mc_logradixl (n, x); }
#	elif MC_TARGET_C11
#	define mcmath_logradix(n, x) _Generic(x \
	, float       : mc_logradixf            \
	, double      : mc_logradix             \
	, long double : mc_logradixl            \
) (mc_cast_expr(unsigned int, n), x)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_logradix(n, x) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_logradixf (mc_cast_expr(unsigned int, n), mc_cast_expr(const float, x))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_logradix  (mc_cast_expr(unsigned int, n), mc_cast_expr(const double, x))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_logradixl (mc_cast_expr(unsigned int, n), mc_cast_expr(const long double, x)) \
		: 0 \
	))
#	else
#	define mcmath_logradix(n, x) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_logradixf (mc_cast_expr(unsigned int, n), mc_cast_expr(const float, x))       \
		: sizeof(x) == sizeof(double)      ? mc_logradix  (mc_cast_expr(unsigned int, n), mc_cast_expr(const double, x))      \
		: sizeof(x) == sizeof(long double) ? mc_logradixl (mc_cast_expr(unsigned int, n), mc_cast_expr(const long double, x)) \
		: 0 \
	)
#	endif
#	endif

#	ifndef mcmath_logbase
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_logbase              (const T& x, const int b)           { mc_unused(x); mc_unused(b); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_logbase<float>       (const float& x, const int b)       { return mc_logbasef (x, b);            }
template <>        MC_TARGET_INLINE double      mcmath_logbase<double>      (const double& x, const int b)      { return mc_logbase  (x, b);            }
template <>        MC_TARGET_INLINE long double mcmath_logbase<long double> (const long double& x, const int b) { return mc_logbasel (x, b);            }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_logbase (const float x, const int b)       { return mc_logbasef (x, b); }
MC_TARGET_ALIAS double      mcmath_logbase (const double x, const int b)      { return mc_logbase  (x, b); }
MC_TARGET_ALIAS long double mcmath_logbase (const long double x, const int b) { return mc_logbasel (x, b); }
#	elif MC_TARGET_C11
#	define mcmath_logbase(x, b) _Generic(x \
	, float       : mc_logbasef            \
	, double      : mc_logbase             \
	, long double : mc_logbasel            \
) (x, mc_cast_expr(const int, b))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_logbase(x, b) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_logbasef (mc_cast_expr(const float, x), mc_cast_expr(const int, b))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_logbase  (mc_cast_expr(const double, x), mc_cast_expr(const int, b))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_logbasel (mc_cast_expr(const long double, x), mc_cast_expr(const int, b)) \
		: 0 \
	))
#	else
#	define mcmath_logbase(x, b) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_logbasef (mc_cast_expr(const float, x), mc_cast_expr(const int, b))       \
		: sizeof(x) == sizeof(double)      ? mc_logbase  (mc_cast_expr(const double, x), mc_cast_expr(const int, b))      \
		: sizeof(x) == sizeof(long double) ? mc_logbasel (mc_cast_expr(const long double, x), mc_cast_expr(const int, b)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_lerp -

#	ifndef mcmath_lerp
#	if MC_TARGET_CPP20 && MC_TARGET_HAVE_LERPFN
#	define mcmath_lerpf ::std::lerp
#	define mcmath_lerp  ::std::lerp
#	define mcmath_lerpl ::std::lerp
#	else
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mcmath_lerp              (const T& x, const T& y, const T& z)                               { mc_unused(x); mc_unused(y); mc_unused(z); return 0; }
template <>        MC_TARGET_INLINE float       mcmath_lerp<float>       (const float& x, const float& y, const float& z)                   { return mc_lerpf (x, y, z);                          }
template <>        MC_TARGET_INLINE double      mcmath_lerp<double>      (const double& x, const double& y, const double& z)                { return mc_lerp  (x, y, z);                          }
template <>        MC_TARGET_INLINE long double mcmath_lerp<long double> (const long double& x, const long double& y, const long double& z) { return mc_lerpl (x, y, z);                          }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mcmath_lerp (const float x, const float y, const float z)                   { return mc_lerpf (x, y, z); }
MC_TARGET_ALIAS double      mcmath_lerp (const double x, const double y, const double z)                { return mc_lerp  (x, y, z); }
MC_TARGET_ALIAS long double mcmath_lerp (const long double x, const long double y, const long double z) { return mc_lerpl( x, y, z); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_lerp(x, y, z) _Generic(x \
	, float       : mc_lerpf               \
	, double      : mc_lerp                \
	, long double : mc_lerpl               \
) (x, mc_cast_expr(MC_TARGET_TYPEOF(x), y), mc_cast_expr(MC_TARGET_TYPEOF(x), z))
#	elif MC_TARGET_C11
#	define mcmath_lerp(x, y, z) _Generic((x)+(y)+(z) \
	, float       : mc_lerpf                         \
	, double      : mc_lerp                          \
	, long double : mc_lerpl                         \
) ((x), (y), (z))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_lerp(x, y, z) mc_cast(MC_TARGET_TYPEOF(x), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? mc_lerpf (mc_cast_expr(const float, x), mc_cast_expr(const float, y), mc_cast_expr(const float, z))                   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? mc_lerp  (mc_cast_expr(const double, x), mc_cast_expr(const double, y), mc_cast_expr(const double, z))                \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? mc_lerpl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y), mc_cast_expr(const long double, z)) \
		: 0 \
	))
#	else
#	define mcmath_lerp(x, y, z) \
	( \
		  sizeof(x) == sizeof(float)       ? mc_lerpf (mc_cast_expr(const float, x), mc_cast_expr(const float, y), mc_cast_expr(const float, z))                   \
		: sizeof(x) == sizeof(double)      ? mc_lerp  (mc_cast_expr(const double, x), mc_cast_expr(const double, y), mc_cast_expr(const double, z))                \
		: sizeof(x) == sizeof(long double) ? mc_lerpl (mc_cast_expr(const long double, x), mc_cast_expr(const long double, y), mc_cast_expr(const long double, z)) \
		: 0 \
	)
#	endif
#	endif
#	endif

#pragma mark - mcmath_cmul -

#	ifndef mcmath_cmul
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T                        mcmath_cmul                           (const T& a, const T& b)                                                 { mc_unused(a); mc_unused(b); return 0; }
template <>        MC_TARGET_INLINE mc_complex_float_t       mcmath_cmul<mc_complex_float_t>       (const mc_complex_float_t & a, const mc_complex_float_t & b)             { return mc_cmulf (a, b);               }
template <>        MC_TARGET_INLINE mc_complex_double_t      mcmath_cmul<mc_complex_double_t>      (const mc_complex_double_t & a, const mc_complex_double_t & b)           { return mc_cmul  (a, b);               }
template <>        MC_TARGET_INLINE mc_complex_long_double_t mcmath_cmul<mc_complex_long_double_t> (const mc_complex_long_double_t & a, const mc_complex_long_double_t & b) { return mc_cmull (a, b);               }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS mc_complex_float_t       mcmath_cmul (const mc_complex_float_t a, const mc_complex_float_t b)             { return mc_cmulf (a, b); }
MC_TARGET_ALIAS mc_complex_double_t      mcmath_cmul (const mc_complex_double_t a, const mc_complex_double_t b)           { return mc_cmul  (a, b); }
MC_TARGET_ALIAS mc_complex_long_double_t mcmath_cmul (const mc_complex_long_double_t a, const mc_complex_long_double_t b) { return mc_cmull (a, b); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_cmul(a, b) _Generic(a   \
	, mc_complex_float_t       : mc_cmulf \
	, mc_complex_double_t      : mc_cmul  \
	, mc_complex_long_double_t : mc_cmull \
) (a, mc_cast_expr(MC_TARGET_TYPEOF(a), b))
#	elif MC_TARGET_C11
#	define mcmath_cmul(a, b) _Generic((a)+(b) \
	, mc_complex_float_t       : mc_cmulf     \
	, mc_complex_double_t      : mc_cmul      \
	, mc_complex_long_double_t : mc_cmull     \
) ((a), (b))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_cmul(a, b) mc_cast(MC_TARGET_TYPEOF(a), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_float_t)       ? mc_cmulf (mc_cast_expr(const mc_complex_float_t, a), mc_cast_expr(const mc_complex_float_t, b))              \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_double_t)      ? mc_cmul  (mc_cast_expr(const mc_complex_double_t, a), mc_cast_expr(const mc_complex_double_t, b))            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_long_double_t) ? mc_cmull (mc_cast_expr(const mc_complex_long_double_t, a), mc_cast_expr(const lmc_complex_long_double_t, b)) \
		: 0 \
	))
#	else
#	define mcmath_cmul(a, b) \
	( \
		  sizeof(x) == sizeof(mc_complex_float_t)       ? mc_cmulf (mc_cast_expr(const mc_complex_float_t, a), mc_cast_expr(const mc_complex_float_t, b))              \
		: sizeof(x) == sizeof(mc_complex_double_t)      ? mc_cmul  (mc_cast_expr(const mc_complex_double_t, a), mc_cast_expr(const mc_complex_double_t, b))            \
		: sizeof(x) == sizeof(mc_complex_long_double_t) ? mc_cmull (mc_cast_expr(const mc_complex_long_double_t, a), mc_cast_expr(const lmc_complex_long_double_t, b)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_cdiv -

#	ifndef mcmath_cdiv
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T                        mcmath_cdiv                           (const T& a, const T& b)                                                 { mc_unused(a); mc_unused(b); return 0; }
template <>        MC_TARGET_INLINE mc_complex_float_t       mcmath_cdiv<mc_complex_float_t>       (const mc_complex_float_t & a, const mc_complex_float_t & b)             { return mc_cdivf (a, b);               }
template <>        MC_TARGET_INLINE mc_complex_double_t      mcmath_cdiv<mc_complex_double_t>      (const mc_complex_double_t & a, const mc_complex_double_t & b)           { return mc_cdiv  (a, b);               }
template <>        MC_TARGET_INLINE mc_complex_long_double_t mcmath_cdiv<mc_complex_long_double_t> (const mc_complex_long_double_t & a, const mc_complex_long_double_t & b) { return mc_cdivl (a, b);               }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS mc_complex_float_t       mcmath_cdiv (const mc_complex_float_t a, const mc_complex_float_t b)             { return mc_cdivf (a, b); }
MC_TARGET_ALIAS mc_complex_double_t      mcmath_cdiv (const mc_complex_double_t a, const mc_complex_double_t b)           { return mc_cdiv  (a, b); }
MC_TARGET_ALIAS mc_complex_long_double_t mcmath_cdiv (const mc_complex_long_double_t a, const mc_complex_long_double_t b) { return mc_cdivl (a, b); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_cdiv(a, b) _Generic(a   \
	, mc_complex_float_t       : mc_cdivf \
	, mc_complex_double_t      : mc_cdiv  \
	, mc_complex_long_double_t : mc_cdivl \
) (a, mc_cast_expr(MC_TARGET_TYPEOF(a), b))
#	elif MC_TARGET_C11
#	define mcmath_cdiv(a, b) _Generic((a)+(b) \
	, mc_complex_float_t       : mc_cdivf     \
	, mc_complex_double_t      : mc_cdiv      \
	, mc_complex_long_double_t : mc_cdivl     \
) ((a), (b))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_cdiv(a, b) mc_cast(MC_TARGET_TYPEOF(a), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_float_t)       ? mc_cdivf (mc_cast_expr(const mc_complex_float_t, a), mc_cast_expr(const mc_complex_float_t, b))              \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_double_t)      ? mc_cdiv  (mc_cast_expr(const mc_complex_double_t, a), mc_cast_expr(const mc_complex_double_t, b))            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_long_double_t) ? mc_cdivl (mc_cast_expr(const mc_complex_long_double_t, a), mc_cast_expr(const lmc_complex_long_double_t, b)) \
		: 0 \
	))
#	else
#	define mcmath_cdiv(a, b) \
	( \
		  sizeof(x) == sizeof(mc_complex_float_t)       ? mc_cdivf (mc_cast_expr(const mc_complex_float_t, a), mc_cast_expr(const mc_complex_float_t, b))              \
		: sizeof(x) == sizeof(mc_complex_double_t)      ? mc_cdiv  (mc_cast_expr(const mc_complex_double_t, a), mc_cast_expr(const mc_complex_double_t, b))            \
		: sizeof(x) == sizeof(mc_complex_long_double_t) ? mc_cdivl (mc_cast_expr(const mc_complex_long_double_t, a), mc_cast_expr(const lmc_complex_long_double_t, b)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_cadd -

#	ifndef mcmath_cadd
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T                        mcmath_cadd                           (const T& a, const T& b)                                                 { mc_unused(a); mc_unused(b); return 0; }
template <>        MC_TARGET_INLINE mc_complex_float_t       mcmath_cadd<mc_complex_float_t>       (const mc_complex_float_t & a, const mc_complex_float_t & b)             { return mc_caddf (a, b);               }
template <>        MC_TARGET_INLINE mc_complex_double_t      mcmath_cadd<mc_complex_double_t>      (const mc_complex_double_t & a, const mc_complex_double_t & b)           { return mc_cadd  (a, b);               }
template <>        MC_TARGET_INLINE mc_complex_long_double_t mcmath_cadd<mc_complex_long_double_t> (const mc_complex_long_double_t & a, const mc_complex_long_double_t & b) { return mc_caddl (a, b);               }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS mc_complex_float_t       mcmath_cadd (const mc_complex_float_t a, const mc_complex_float_t b)             { return mc_caddf (a, b); }
MC_TARGET_ALIAS mc_complex_double_t      mcmath_cadd (const mc_complex_double_t a, const mc_complex_double_t b)           { return mc_cadd  (a, b); }
MC_TARGET_ALIAS mc_complex_long_double_t mcmath_cadd (const mc_complex_long_double_t a, const mc_complex_long_double_t b) { return mc_caddl (a, b); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_cadd(a, b) _Generic(a   \
	, mc_complex_float_t       : mc_caddf \
	, mc_complex_double_t      : mc_cadd  \
	, mc_complex_long_double_t : mc_caddl \
) (a, mc_cast_expr(MC_TARGET_TYPEOF(a), b))
#	elif MC_TARGET_C11
#	define mcmath_cadd(a, b) _Generic((a)+(b) \
	, mc_complex_float_t       : mc_caddf     \
	, mc_complex_double_t      : mc_cadd      \
	, mc_complex_long_double_t : mc_caddl     \
) ((a), (b))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_cadd(a, b) mc_cast(MC_TARGET_TYPEOF(a), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_float_t)       ? mc_caddf (mc_cast_expr(const mc_complex_float_t, a), mc_cast_expr(const mc_complex_float_t, b))              \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_double_t)      ? mc_cadd  (mc_cast_expr(const mc_complex_double_t, a), mc_cast_expr(const mc_complex_double_t, b))            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_long_double_t) ? mc_caddl (mc_cast_expr(const mc_complex_long_double_t, a), mc_cast_expr(const lmc_complex_long_double_t, b)) \
		: 0 \
	))
#	else
#	define mcmath_cadd(a, b) \
	( \
		  sizeof(x) == sizeof(mc_complex_float_t)       ? mc_caddf (mc_cast_expr(const mc_complex_float_t, a), mc_cast_expr(const mc_complex_float_t, b))              \
		: sizeof(x) == sizeof(mc_complex_double_t)      ? mc_cadd  (mc_cast_expr(const mc_complex_double_t, a), mc_cast_expr(const mc_complex_double_t, b))            \
		: sizeof(x) == sizeof(mc_complex_long_double_t) ? mc_caddl (mc_cast_expr(const mc_complex_long_double_t, a), mc_cast_expr(const lmc_complex_long_double_t, b)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_csub -

#	ifndef mcmath_csub
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T                        mcmath_csub                           (const T& a, const T& b)                                                 { mc_unused(a); mc_unused(b); return 0; }
template <>        MC_TARGET_INLINE mc_complex_float_t       mcmath_csub<mc_complex_float_t>       (const mc_complex_float_t & a, const mc_complex_float_t & b)             { return mc_csubf (a, b);               }
template <>        MC_TARGET_INLINE mc_complex_double_t      mcmath_csub<mc_complex_double_t>      (const mc_complex_double_t & a, const mc_complex_double_t & b)           { return mc_csub  (a, b);               }
template <>        MC_TARGET_INLINE mc_complex_long_double_t mcmath_csub<mc_complex_long_double_t> (const mc_complex_long_double_t & a, const mc_complex_long_double_t & b) { return mc_csubl (a, b);               }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS mc_complex_float_t       mcmath_csub (const mc_complex_float_t a, const mc_complex_float_t b)             { return mc_csubf (a, b); }
MC_TARGET_ALIAS mc_complex_double_t      mcmath_csub (const mc_complex_double_t a, const mc_complex_double_t b)           { return mc_csub  (a, b); }
MC_TARGET_ALIAS mc_complex_long_double_t mcmath_csub (const mc_complex_long_double_t a, const mc_complex_long_double_t b) { return mc_csubl (a, b); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_csub(a, b) _Generic(a   \
	, mc_complex_float_t       : mc_csubf \
	, mc_complex_double_t      : mc_csub  \
	, mc_complex_long_double_t : mc_csubl \
) (a, mc_cast_expr(MC_TARGET_TYPEOF(a), b))
#	elif MC_TARGET_C11
#	define mcmath_csub(a, b) _Generic((a)+(b) \
	, mc_complex_float_t       : mc_csubf     \
	, mc_complex_double_t      : mc_csub      \
	, mc_complex_long_double_t : mc_csubl     \
) ((a), (b))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_csub(a, b) mc_cast(MC_TARGET_TYPEOF(a), \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_float_t)       ? mc_csubf (mc_cast_expr(const mc_complex_float_t, a), mc_cast_expr(const mc_complex_float_t, b))              \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_double_t)      ? mc_csub  (mc_cast_expr(const mc_complex_double_t, a), mc_cast_expr(const mc_complex_double_t, b))            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_long_double_t) ? mc_csubl (mc_cast_expr(const mc_complex_long_double_t, a), mc_cast_expr(const lmc_complex_long_double_t, b)) \
		: 0 \
	))
#	else
#	define mcmath_csub(a, b) \
	( \
		  sizeof(x) == sizeof(mc_complex_float_t)       ? mc_csubf (mc_cast_expr(const mc_complex_float_t, a), mc_cast_expr(const mc_complex_float_t, b))              \
		: sizeof(x) == sizeof(mc_complex_double_t)      ? mc_csub  (mc_cast_expr(const mc_complex_double_t, a), mc_cast_expr(const mc_complex_double_t, b))            \
		: sizeof(x) == sizeof(mc_complex_long_double_t) ? mc_csubl (mc_cast_expr(const mc_complex_long_double_t, a), mc_cast_expr(const lmc_complex_long_double_t, b)) \
		: 0 \
	)
#	endif
#	endif

#pragma mark - mcmath_ciseq -

#	ifndef mcmath_ciseq
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE int mcmath_ciseq                           (const T& a, const T& b)                                                 { mc_unused(a); mc_unused(b); return 0; }
template <>        MC_TARGET_INLINE int mcmath_ciseq<mc_complex_float_t>       (const mc_complex_float_t & a, const mc_complex_float_t & b)             { return mc_ciseqf (a, b);              }
template <>        MC_TARGET_INLINE int mcmath_ciseq<mc_complex_double_t>      (const mc_complex_double_t & a, const mc_complex_double_t & b)           { return mc_ciseq  (a, b);              }
template <>        MC_TARGET_INLINE int mcmath_ciseq<mc_complex_long_double_t> (const mc_complex_long_double_t & a, const mc_complex_long_double_t & b) { return mc_ciseql (a, b);              }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS int mcmath_ciseq (const mc_complex_float_t a, const mc_complex_float_t b)             { return mc_ciseqf (a, b); }
MC_TARGET_ALIAS int mcmath_ciseq (const mc_complex_double_t a, const mc_complex_double_t b)           { return mc_ciseq  (a, b); }
MC_TARGET_ALIAS int mcmath_ciseq (const mc_complex_long_double_t a, const mc_complex_long_double_t b) { return mc_ciseql (a, b); }
#	elif MC_TARGET_C11 && MC_TARGET_HAVE_TYPEOF
#	define mcmath_ciseq(a, b) _Generic(a   \
	, mc_complex_float_t       : mc_ciseqf \
	, mc_complex_double_t      : mc_ciseq  \
	, mc_complex_long_double_t : mc_ciseql \
) (a, mc_cast_expr(MC_TARGET_TYPEOF(a), b))
#	elif MC_TARGET_C11
#	define mcmath_ciseq(a, b) _Generic((a)+(b) \
	, mc_complex_float_t       : mc_ciseqf     \
	, mc_complex_double_t      : mc_ciseq      \
	, mc_complex_long_double_t : mc_ciseql     \
) ((a), (b))
#	elif MC_TARGET_HAVE_TYPEOF
#	define mcmath_ciseq(a, b) \
	( \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_float_t)       ? mc_ciseqf (mc_cast_expr(const mc_complex_float_t, a), mc_cast_expr(const mc_complex_float_t, b))              \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_double_t)      ? mc_ciseq  (mc_cast_expr(const mc_complex_double_t, a), mc_cast_expr(const mc_complex_double_t, b))            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(a), mc_complex_long_double_t) ? mc_ciseql (mc_cast_expr(const mc_complex_long_double_t, a), mc_cast_expr(const lmc_complex_long_double_t, b)) \
		: 0 \
	)
#	else
#	define mcmath_ciseq(a, b) \
	( \
		  sizeof(x) == sizeof(mc_complex_float_t)       ? mc_ciseqf (mc_cast_expr(const mc_complex_float_t, a), mc_cast_expr(const mc_complex_float_t, b))              \
		: sizeof(x) == sizeof(mc_complex_double_t)      ? mc_ciseq  (mc_cast_expr(const mc_complex_double_t, a), mc_cast_expr(const mc_complex_double_t, b))            \
		: sizeof(x) == sizeof(mc_complex_long_double_t) ? mc_ciseql (mc_cast_expr(const mc_complex_long_double_t, a), mc_cast_expr(const lmc_complex_long_double_t, b)) \
		: 0 \
	)
#	endif
#	endif

#endif /* !MCMATH_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_target.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#ifndef MC_TARGET_H
#define MC_TARGET_H

#	if    (defined(__INTEL_COMPILER))                                    \
		|| (defined(__ICL))                                               \
		|| (defined(__ibmxl__))                                           \
		|| (defined(__xlc__))                                             \
		|| (defined(__xlC__))                                             \
		|| (defined(_MSC_VER)         && (defined(__STDC__) && __STDC__)) \
		|| (defined(_MSC_VER)         && _MSC_VER < 1916)
#		error "C99 or CPP11 compiler and Posix 1-2001 CRT required."
#	else
#		undef  MC_DISABLE_TYPEOF
#		undef  MC_DISABLE_TGMATH
#		undef  MC_DISABLE_ALLOW_CPP_CMATH
#		define MC_DISABLE_TYPEOF 1
#		define MC_DISABLE_TGMATH 1
#		define MC_DISABLE_ALLOW_CPP_CMATH 1
#	endif

#	if defined(__GLIBC__)
#		ifndef _GNU_SOURCE
#			define _GNU_SOURCE
#		endif
#	endif

#	if  defined(__unix__)      \
	||  defined(__linux__)     \
	||  defined(__gnu_linux__) \
	||  defined(__bsdi__)      \
	||  defined(__FreeBSD__)   \
	||  defined(__NetBSD__)    \
	||  defined(__OpenBSD__)   \
	||  defined(__DragonFly__) \
	|| (defined(__APPLE__) && defined(__MACH__))
#	include <unistd.h>
#	endif

#	if defined(__APPLE__) && defined(__MACH__)
#		ifndef __MATH_LONG_DOUBLE_CONSTANTS
#			define __MATH_LONG_DOUBLE_CONSTANTS 1
#		endif
#	endif

#	if MC_DISABLE_REENTRANT
#		undef  _REENTRANT
#	else
#		undef  _REENTRANT
#		define _REENTRANT 1
#	endif

#	if MC_DISABLE_TGMATH
#		undef  MC_TARGET_HAVE_TGMATH
#		define MC_TARGET_HAVE_TGMATH 0
#	else
#		undef  MC_TARGET_HAVE_TGMATH
#		define MC_TARGET_HAVE_TGMATH 1
#	endif

#	if MC_DISABLE_ALLOW_CPP_CMATH
#		undef  MC_TARGET_ALLOW_CPP_CMATH
#		define MC_TARGET_ALLOW_CPP_CMATH 0
#	else
#		undef  MC_TARGET_ALLOW_CPP_CMATH
#		define MC_TARGET_ALLOW_CPP_CMATH 1
#	endif

#	if MC_DISABLE_LOG2
#		undef  MC_TARGET_HAVE_LOG2
#		define MC_TARGET_HAVE_LOG2 0
#	else
#		undef  MC_TARGET_HAVE_LOG2
#		define MC_TARGET_HAVE_LOG2 1
#	endif

#	if MC_DISABLE_INLINE
#		undef  MC_TARGET_INLINE
#		define MC_TARGET_INLINE inline
#	endif

#	if defined(__STDC__) && !defined(__cplusplus)
#		define MC_TARGET_C89 1
#		if defined(__STDC_VERSION__)
#			define MC_TARGET_C90 1
#			if (__STDC_VERSION__ >= 199409L)
#				define MC_TARGET_C94 1
#			endif
#			if (__STDC_VERSION__ >= 199901L)
#				define MC_TARGET_C99 1
#			endif
#			if (__STDC_VERSION__ >= 201112L)
#				define MC_TARGET_C11 1
#			endif
#			if (__STDC_VERSION__ >= 201710L)
#				define MC_TARGET_C17 1
#			endif
#		endif
#	endif

#	if defined(__cplusplus)
#		undef MC_TARGET_C89
#		undef MC_TARGET_C90
#		undef MC_TARGET_C99
#		undef MC_TARGET_C11
#		undef MC_TARGET_C17
#		define MC_TARGET_CPP98 1
#		if (__cplusplus >= 201103L)
#			define MC_TARGET_CPP11 1
#		endif
#		if (__cplusplus >= 201402L)
#			define MC_TARGET_CPP14 1
#		endif
#		if (__cplusplus >= 201703L)
#			define MC_TARGET_CPP17 1
#		endif
#		if (__cplusplus > 201703L)
#			define MC_TARGET_CPP20 1
#		endif
#	endif

#	undef MC_TARGET_MSVC_CPP
#	if defined(_MSC_VER)
#		if !defined(__cplusplus)
#			error "MSC C++ support only."
#		endif
#		define MC_TARGET_MSVC_CPP 1
#		ifndef MC_TARGET_CPP98
#			define MC_TARGET_CPP98 1
#		endif
#		ifndef MC_TARGET_CPP11
#			define MC_TARGET_CPP11 1
#		endif
#		ifndef MC_TARGET_CPP14
#			define MC_TARGET_CPP14 1
#		endif
#	endif

#	if MC_TARGET_MSVC_CPP
#		undef  _USE_MATH_DEFINES
#		define _USE_MATH_DEFINES 1
#	endif

#	if MC_TARGET_CPP98 && !MC_TARGET_CPP11
#		ifndef __STDC_CONSTANT_MACROS
#			define __STDC_CONSTANT_MACROS
#		endif
#		ifndef __STDC_LIMIT_MACROS
#			define __STDC_LIMIT_MACROS
#		endif
#		ifndef __STDC_FORMAT_MACROS
#			define __STDC_FORMAT_MACROS
#		endif
#	endif

#	if MC_TARGET_C89 && !MC_TARGET_C99
#	if defined(__GNUC__)
#		define inline __inline__
#	else
#		define inline __inline
#	endif
#	endif

#	if MC_TARGET_C99
#		define MC_TARGET_RESTRICT restrict
#	else
#		define MC_TARGET_RESTRICT
#	endif

#	if !defined(MC_TARGET_INLINE)
#	if MC_TARGET_C99 || MC_TARGET_CPP98
#		if (((defined(__GNUC__) && __GNUC__ + 0 >= 4)) || defined(__clang__))
#			define MC_TARGET_INLINE __inline__ __attribute__((__always_inline__, __unused__))
#		elif defined(__GNUC__)
#			define MC_TARGET_INLINE __inline__
#		elif defined(_MSC_VER)
#			define MC_TARGET_INLINE __forceinline
#		else
#			define MC_TARGET_INLINE inline
#		endif
#	else
#		define MC_TARGET_INLINE
#	endif
#	endif

#	if !defined(MC_TARGET_PROC)
#		define MC_TARGET_PROC static MC_TARGET_INLINE
#	endif

#	if !defined(MC_TARGET_FUNC)
#		define MC_TARGET_FUNC static MC_TARGET_INLINE
#	endif

#	if MC_DISABLE_OVERLOADABLE
#		undef  MC_TARGET_HAVE_OVERLOADABLE
#		undef  MC_TARGET_OVERLOADABLE
#		undef  MC_TARGET_ALIAS
#		define MC_TARGET_HAVE_OVERLOADABLE   0
#		define MC_TARGET_OVERLOADABLE        MC_TARGET_INLINE
#		define MC_TARGET_ALIAS               static MC_TARGET_OVERLOADABLE
#	else
#	if defined(__clang__) && __has_attribute(__overloadable__)
#		undef  MC_TARGET_HAVE_OVERLOADABLE
#		undef  MC_TARGET_OVERLOADABLE
#		undef  MC_TARGET_ALIAS
#		define MC_TARGET_HAVE_OVERLOADABLE   1
#		define MC_TARGET_OVERLOADABLE        __attribute__((__overloadable__, __always_inline__, __unused__))
#		define MC_TARGET_ALIAS               static MC_TARGET_OVERLOADABLE
#	else
#		undef  MC_TARGET_HAVE_OVERLOADABLE
#		undef  MC_TARGET_OVERLOADABLE
#		undef  MC_TARGET_ALIAS
#		define MC_TARGET_HAVE_OVERLOADABLE   0
#		define MC_TARGET_OVERLOADABLE
#		define MC_TARGET_ALIAS
#	endif
#	endif

#	if defined(WIN32) && (defined(_MSC_VER) || defined(__ICL))
#		define MC_TARGET_THREAD_LOCAL __declspec(thread)
#	elif defined(__clang__)
#		if __has_feature(c_thread_local) || __has_extension(c_thread_local)
#			define MC_TARGET_THREAD_LOCAL _Thread_local
#		else
#			define MC_TARGET_THREAD_LOCAL
#		endif
#	elif defined(__GNUG__) && (__GNUC__ <= 4 && __GNUC_MINOR__ < 80)
#		define MC_TARGET_THREAD_LOCAL __thread
#	elif defined(__GNUG__) && (__GNUC__ == 4 && __GNUC_MINOR__ >= 80)
#		define MC_TARGET_THREAD_LOCAL _Thread_local
#	elif defined(__GNUG__) && (__GNUC__ + 0 > 4)
#		if MC_TARGET_CPP11
#			define MC_TARGET_THREAD_LOCAL thread_local
#		else
#			define MC_TARGET_THREAD_LOCAL _Thread_local
#		endif
#	elif (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L) && !defined(__STDC_NO_THREADS__)
#		define MC_TARGET_THREAD_LOCAL _Thread_local
#	else
#		define MC_TARGET_THREAD_LOCAL
#	endif

#	if MC_DISABLE_TYPEOF
#		undef  MC_TARGET_HAVE_TYPEOF
#		undef  MC_TARGET_HAVE_AUTOTYPE
#		undef  MC_TARGET_AUTOTYPE
#		undef  MC_TARGET_TYPEOF
#		undef  MC_TARGET_TYPEISOF
#		define MC_TARGET_HAVE_TYPEOF      0
#		define MC_TARGET_HAVE_AUTOTYPE    0
#		define MC_TARGET_AUTOTYPE         void *
#		define MC_TARGET_TYPEOF(X)        void *
#		define MC_TARGET_TYPEISOF(T1, T2) 0
#	else
#	if MC_TARGET_C89 || MC_TARGET_CPP98
#		undef  MC_TARGET_HAVE_AUTOTYPE
#		undef  MC_TARGET_HAVE_TYPEOF
#		undef  MC_TARGET_AUTOTYPE
#		undef  MC_TARGET_TYPEOF
#		undef  MC_TARGET_TYPEISOF
#		define MC_TARGET_HAVE_TYPEOF 1
#		if MC_TARGET_CPP11
#			define MC_TARGET_HAVE_AUTOTYPE 1
#			define MC_TARGET_AUTOTYPE      auto
#		elif ((defined(__GNUC__) && (__GNUC__ + 0 >= 5)) || defined(__clang__))
#			define MC_TARGET_HAVE_AUTOTYPE 1
#			define MC_TARGET_AUTOTYPE      __auto_type
#		endif
#		if ((defined(__GNUC__) && (__GNUC__ + 0 >= 3)) || defined(__clang__))
#			define MC_TARGET_TYPEISOF(type1, type2) ((__builtin_types_compatible_p(type1, type2)) ? 1 : 0)
#		elif MC_TARGET_CPP11
#			define MC_TARGET_TYPEISOF(type1, type2) \
			( \
				::std::is_same< \
					  typename ::std::decay<type1>::type \
					, typename ::std::decay<type2>::type \
				>::value == true ? 1 : 0 \
			)
#		elif MC_TARGET_CPP98
			template<class T>        struct mc_remove_const                { typedef T type; };
			template<class T>        struct mc_remove_const<const T>       { typedef T type; };
			template<class T>        struct mc_remove_volatile             { typedef T type; };
			template<class T>        struct mc_remove_volatile<volatile T> { typedef T type; };
			template<class T>        struct mc_remove_cv
			{
				typedef typename mc_remove_volatile<typename mc_remove_const<T>::type>::type type;
			};
			template<class T>        struct mc_remove_reference            { typedef T type; };
			template<class T>        struct mc_remove_reference<T&>        { typedef T type; };
			template<class T>        struct mc_decay
			{
				typedef typename mc_remove_reference<typename mc_remove_cv<T>::type>::type type;
			};
			template<class T, class U> struct mc_is_same                   { static const bool value = false; };
			template<class T>          struct mc_is_same<T, T>             { static const bool value = true;  };
#			define MC_TARGET_TYPEISOF(type1, type2) \
			( \
				mc_is_same< \
					  typename mc_decay<type1>::type \
					, typename mc_decay<type2>::type \
				>::value == true ? 1 : 0 \
			)
#		endif
#		if MC_TARGET_MSVC_CPP
#			define MC_TARGET_TYPEOF(x) ::std::decay<decltype((x))>::type
#		else
#			define MC_TARGET_TYPEOF(x) __typeof__((x) + 0)
#		endif
#	else
#		undef  MC_TARGET_HAVE_TYPEOF
#		undef  MC_TARGET_HAVE_AUTOTYPE
#		undef  MC_TARGET_AUTOTYPE
#		undef  MC_TARGET_TYPEOF
#		undef  MC_TARGET_TYPEISOF
#		define MC_TARGET_HAVE_TYPEOF      0
#		define MC_TARGET_HAVE_AUTOTYPE    0
#		define MC_TARGET_AUTOTYPE         void *
#		define MC_TARGET_TYPEOF(X)        void *
#		define MC_TARGET_TYPEISOF(T1, T2) 0
#	endif
#	endif

#	if MC_TARGET_CPP98
#		include <cassert>
#	else
#		include <assert.h>
#	endif

#	if MC_TARGET_CPP11
#	define mc_static_assert(COND, IDF) ::static_assert(COND, #IDF)
#	elif MC_TARGET_C11
#	define mc_static_assert(COND, IDF) _Static_assert(COND, #IDF)
#	else
#	define mc_static_assert(COND, IDF) typedef char static_assertion_##IDF[(COND) ? 1 : -1]
#	endif

#	if MC_TARGET_CPP98
#	define mc_cast(t, x)      static_cast<t>(x)
#	define mc_cast_expr(t, x) static_cast<t>((x))
#	if MC_TARGET_CPP11
#		define MC_NULLPTR nullptr
#	else
#		define MC_NULLPTR NULL
#	endif
#	else
#	define mc_cast(t, x) (t)x
#	define mc_cast_expr(t, x) mc_cast(t, (x))
#	define MC_NULLPTR NULL
#	endif

#	define mc_nonnullptr(p) ((p) != MC_NULLPTR)
#	define mc_unused(x)     (void)x

#	define mc_scope_begin do     {
#	define mc_scope_end   break; } while (0)

#	undef MC_TARGET_LONG_IS64
#	if defined(INT_TYPE_SIZE) && defined(LONG_TYPE_SIZE) && defined(LONG_LONG_TYPE_SIZE)
#		if (INT_TYPE_SIZE < LONG_TYPE_SIZE && !(LONG_TYPE_SIZE < LONG_LONG_TYPE_SIZE)
#			define MC_TARGET_LONG_64BIT 1
#		endif
#	endif

#	if MC_TARGET_CPP98
#		include <cstdio>
#		include <climits>
#		include <ctime>
#		include <cstdlib>
#		include <cstring>
#		include <cctype>
#		include <cstdint>
#		include <cinttypes>
#		include <cfenv>
#		include <cfloat>
#		include <cmath>
#		include <algorithm>
#	elif MC_TARGET_HAVE_TGMATH
#		include <stdio.h>
#		include <limits.h>
#		include <time.h>
#		include <stdlib.h>
#		include <string.h>
#		include <ctype.h>
#		include <stdint.h>
#		include <inttypes.h>
#		include <fenv.h>
#		include <float.h>
#		include <tgmath.h>
#	else
#		include <stdio.h>
#		include <limits.h>
#		include <time.h>
#		include <stdlib.h>
#		include <string.h>
#		include <ctype.h>
#		include <stdint.h>
#		include <inttypes.h>
#		include <fenv.h>
#	if !__STDC_NO_COMPLEX__
#		include <complex.h>
#	endif
#		include <float.h>
#		include <math.h>
#	endif

#	if defined(__APPLE__)
	extern float       lgammaf_r(float, int *);
	extern double      lgamma_r(double, int *);
	extern long double lgammal_r(long double, int *);
#	if MC_TARGET_CPP98
	namespace std {
		using ::lgammaf_r;
		using ::lgamma_r;
		using ::lgammal_r;
	}
#	endif
#	endif

#	undef  MC_TARGET_HAVE_FMA
#	define MC_TARGET_HAVE_FMA 1

#	if DBL_MANT_DIG < LDBL_MANT_DIG
#		define MC_TARGET_HAVE_LONG_DOUBLE 1
#	else
#		define MC_TARGET_HAVE_LONG_DOUBLE 0
#	endif

#	ifndef FLT_DECIMAL_DIG
#		ifdef __FLT_DECIMAL_DIG__
#			define FLT_DECIMAL_DIG  __FLT_DECIMAL_DIG__
#		else
#			define FLT_DECIMAL_DIG  (FLT_DIG  + 3)
#		endif
#	endif
#	ifndef DBL_DECIMAL_DIG
#		ifdef __DBL_DECIMAL_DIG__
#			define DBL_DECIMAL_DIG  __DBL_DECIMAL_DIG__
#		else
#			define DBL_DECIMAL_DIG  (DBL_DIG  + 2)
#		endif
#	endif
#	ifndef LDBL_DECIMAL_DIG
#		ifdef __LDBL_DECIMAL_DIG__
#			define LDBL_DECIMAL_DIG  __LDBL_DECIMAL_DIG__
#		else
#			if MC_TARGET_HAVE_LONG_DOUBLE
#				define LDBL_DECIMAL_DIG (LDBL_DIG + 3)
#			else
#				define LDBL_DECIMAL_DIG (LDBL_DIG + 2)
#			endif
#		endif
#	endif

#	define MC_TARGET_LONG_DOUBLE_ALIAS 11
#	define MC_TARGET_LONG_DOUBLE_X87   12
#	define MC_TARGET_LONG_DOUBLE_IBM   13
#	define MC_TARGET_LONG_DOUBLE_IEEE  14

#	if   LDBL_MANT_DIG + 0 == 64
#		define MC_TARGET_LONG_DOUBLE_TYPE MC_TARGET_LONG_DOUBLE_X87
#	elif LDBL_MANT_DIG + 0 == 106
#		define MC_TARGET_LONG_DOUBLE_TYPE MC_TARGET_LONG_DOUBLE_IBM
#	elif LDBL_MANT_DIG + 0 == 113
#		define MC_TARGET_LONG_DOUBLE_TYPE MC_TARGET_LONG_DOUBLE_IEEE
#	else
#		define MC_TARGET_LONG_DOUBLE_TYPE MC_TARGET_LONG_DOUBLE_ALIAS
#	endif

#	if MC_TARGET_C99 && !MC_TARGET_BUILTIN_COMPLEX && !__STDC_NO_COMPLEX__
#		if (defined(_Imaginary_I) || defined(_Complex_I)) && !defined(__STDC_IEC_559_COMPLEX__)
#			define __STDC_IEC_559_COMPLEX__ 1
#		endif
#	endif

#	if MC_TARGET_BUILTIN_COMPLEX
#		undef  MC_TARGET_C99_COMPLEX
#		define MC_TARGET_C99_COMPLEX 0
#	endif

#	if MC_TARGET_BLAS_USE_NATIVE

#	undef OPENBLAS_COMPLEX_C99
#	undef OPENBLAS_COMPLEX_STRUCT

#	if MC_TARGET_BLAS_USE_ACCELERATE
#		include <Accelerate/Accelerate.h>
		typedef __CLPK_complex                           mc_complex_float_t;
		typedef __CLPK_doublecomplex                     mc_complex_double_t;
		typedef struct { long double r; long double i; } mc_complex_long_double_t;
#		undef  MC_TARGET_BUILTIN_COMPLEX
#		define MC_TARGET_BUILTIN_COMPLEX 1
#	elif MC_TARGET_BLAS_USE_VECLIB
#		include <vecLib/clapack.h>
#		include <vecLib/cblas.h>
		typedef __CLPK_complex                           mc_complex_float_t;
		typedef __CLPK_doublecomplex                     mc_complex_double_t;
		typedef struct { long double r; long double i; } mc_complex_long_double_t;
#		undef  MC_TARGET_BUILTIN_COMPLEX
#		define MC_TARGET_BUILTIN_COMPLEX 1
#	elif MC_TARGET_BLAS_USE_OPENBLAS
#		if defined __has_include
#			if __has_include("cblas_openblas.h")
#				include "cblas_openblas.h"
#			elif __has_include("cblas-openblas.h")
#				include "cblas-openblas.h"
#			else
#				include "cblas.h"
#			endif
#		else
#			include "cblas.h"
#		endif
#		if defined(OPENBLAS_COMPLEX_C99) || defined(OPENBLAS_COMPLEX_STRUCT)
			typedef openblas_complex_float                            mc_complex_float_t;
			typedef openblas_complex_double                           mc_complex_double_t;
#			if defined(OPENBLAS_COMPLEX_STRUCT)
#				undef  MC_TARGET_BUILTIN_COMPLEX
#				define MC_TARGET_BUILTIN_COMPLEX 1
				typedef struct { long double real; long double imag; } mc_complex_long_double_t;
#			else
				typedef long double _Complex                           mc_complex_long_double_t;
#			endif
#		else
#			undef MC_TARGET_BLAS_USE_NATIVE
#			error "OpenBlas header not found."
#		endif
#	else
#		undef MC_TARGET_BLAS_USE_NATIVE
#		error "Blas native target not found."
#	endif
#	endif

#	if ((!MC_TARGET_CPP98 && !MC_TARGET_BUILTIN_COMPLEX) || defined(OPENBLAS_COMPLEX_C99))
#		if ((MC_TARGET_C99 && defined(__STDC_IEC_559_COMPLEX__) && !__STDC_NO_COMPLEX__) || defined(OPENBLAS_COMPLEX_C99))
#			undef  MC_TARGET_C99_COMPLEX
#			define MC_TARGET_C99_COMPLEX 1
#			if !defined(OPENBLAS_COMPLEX_C99)
#				define  mc_complex(type)        type _Complex
				typedef mc_complex(float)       mc_complex_float_t;
				typedef mc_complex(double)      mc_complex_double_t;
				typedef mc_complex(long double) mc_complex_long_double_t;
#			endif
#			if !MC_TARGET_C11
#			ifndef CMPLXF
#				define CMPLXF(re, im) ((float _Complex)      ((float)(re)       + _Imaginary_I * (float)(im)))
#			endif
#			ifndef CMPLX
#				define CMPLX(re, im)  ((double _Complex)     ((double)(re)      + _Imaginary_I * (double)(im)))
#			endif
#			ifndef CMPLXL
#				define CMPLXL(re, im) ((long double _Complex)((long double)(re) + _Imaginary_I * (long double)(im)))
#			endif
#			endif
#			define mc_cmplxf(re, im) CMPLXF(re, im)
#			define mc_cmplx(re, im)  CMPLX(re, im)
#			define mc_cmplxl(re, im) CMPLXL(re, im)
#			define mc_cmplxrf(c)     crealf(c)
#			define mc_cmplxr(c)      creal (c)
#			define mc_cmplxrl(c)     creall(c)
#			define mc_cmplxif(c)     cimagf(c)
#			define mc_cmplxi(c)      cimag (c)
#			define mc_cmplxil(c)     cimagl(c)
#		endif
#	endif

#	if !MC_TARGET_C99_COMPLEX
#	if MC_TARGET_BLAS_USE_NATIVE
#		if !defined(OPENBLAS_COMPLEX_STRUCT) && !(MC_TARGET_BLAS_USE_ACCELERATE || MC_TARGET_BLAS_USE_VECLIB)
#			define mc_complex(type) struct { type real; type imag; }
#		endif
#	else
#		define mc_complex(type) struct { type u_re; type u_im; }
#	endif
#	if MC_TARGET_BLAS_USE_NATIVE
#		if !defined(OPENBLAS_COMPLEX_STRUCT) && !(MC_TARGET_BLAS_USE_ACCELERATE || MC_TARGET_BLAS_USE_VECLIB)
			typedef mc_complex(float)       mc_complex_float_t;
			typedef mc_complex(double)      mc_complex_double_t;
			typedef mc_complex(long double) mc_complex_long_double_t;
#		endif
#	else
	typedef mc_complex(float)       mc_complex_float_t;
	typedef mc_complex(double)      mc_complex_double_t;
	typedef mc_complex(long double) mc_complex_long_double_t;
#	endif
#	if MC_TARGET_CPP11
#		define mc_cmplxf(re, im) { (float)(re)      , (float)(im)       }
#		define mc_cmplx(re, im)  { (double)(re)     , (double)(im)      }
#		define mc_cmplxl(re, im) { (long double)(re), (long double)(im) }
#	elif MC_TARGET_CPP98
#		if (defined(__GNUC__) || defined(__clang__))
#			define mc_cmplxf(re, im) __extension__ (mc_complex_float_t)       { (float)(re)      , (float)(im)       }
#			define mc_cmplx(re, im)  __extension__ (mc_complex_double_t)      { (double)(re)     , (double)(im)      }
#			define mc_cmplxl(re, im) __extension__ (mc_complex_long_double_t) { (long double)(re), (long double)(im) }
#		else
			MC_TARGET_FUNC mc_complex_float_t       mc_cmplxf (const float re, const float im)              { mc_complex_float_t       z = { re, im }; return z; }
			MC_TARGET_FUNC mc_complex_double_t      mc_cmplx  (const double re, const double im)            { mc_complex_double_t      z = { re, im }; return z; }
			MC_TARGET_FUNC mc_complex_long_double_t mc_cmplxl (const long double re, const long  double im) { mc_complex_long_double_t z = { re, im }; return z; }
#		endif
#	elif MC_TARGET_C99
#		if MC_TARGET_BLAS_USE_NATIVE
#			if MC_TARGET_BLAS_USE_ACCELERATE || MC_TARGET_BLAS_USE_VECLIB
#				define mc_cmplxf(re, im) (mc_complex_float_t)       { .r = (float)(re)      , .i = (float)(im)       }
#				define mc_cmplx(re, im)  (mc_complex_double_t)      { .r = (double)(re)     , .i = (double)(im)      }
#				define mc_cmplxl(re, im) (mc_complex_long_double_t) { .r = (long double)(re), .i = (long double)(im) }
#			else
#				define mc_cmplxf(re, im) (mc_complex_float_t)       { .real = (float)(re)      , .imag = (float)(im)       }
#				define mc_cmplx(re, im)  (mc_complex_double_t)      { .real = (double)(re)     , .imag = (double)(im)      }
#				define mc_cmplxl(re, im) (mc_complex_long_double_t) { .real = (long double)(re), .imag = (long double)(im) }
#			endif
#		else
#			define mc_cmplxf(re, im) (mc_complex_float_t)       { .u_re = (float)(re)      , .u_im = (float)(im)       }
#			define mc_cmplx(re, im)  (mc_complex_double_t)      { .u_re = (double)(re)     , .u_im = (double)(im)      }
#			define mc_cmplxl(re, im) (mc_complex_long_double_t) { .u_re = (long double)(re), .u_im = (long double)(im) }
#		endif
#	elif MC_TARGET_C89
#		define mc_cmplxf(re, im) __extension__ (mc_complex_float_t)       { (float)(re)      , (float)(im)       }
#		define mc_cmplx(re, im)  __extension__ (mc_complex_double_t)      { (double)(re)     , (double)(im)      }
#		define mc_cmplxl(re, im) __extension__ (mc_complex_long_double_t) { (long double)(re), (long double)(im) }
#	endif
#	if MC_TARGET_BLAS_USE_NATIVE
#	if MC_TARGET_BLAS_USE_ACCELERATE || MC_TARGET_BLAS_USE_VECLIB
#		define mc_cmplxrf(c) c.r
#		define mc_cmplxr(c)  c.r
#		define mc_cmplxrl(c) c.r
#		define mc_cmplxif(c) c.i
#		define mc_cmplxi(c)  c.i
#		define mc_cmplxil(c) c.i
#	else
#		define mc_cmplxrf(c) c.real
#		define mc_cmplxr(c)  c.real
#		define mc_cmplxrl(c) c.real
#		define mc_cmplxif(c) c.imag
#		define mc_cmplxi(c)  c.imag
#		define mc_cmplxil(c) c.imag
#	endif
#	else
#		define mc_cmplxrf(c) c.u_re
#		define mc_cmplxr(c)  c.u_re
#		define mc_cmplxrl(c) c.u_re
#		define mc_cmplxif(c) c.u_im
#		define mc_cmplxi(c)  c.u_im
#		define mc_cmplxil(c) c.u_im
#	endif
#	endif

#	if defined(__SSE__) && __SSE__
#		undef  MC_TARGET_HAVE_SSE
#		define MC_TARGET_HAVE_SSE 1
#		include <immintrin.h>
#	endif

#	if (defined(__ARM_NEON__) && __ARM_NEON__) || (defined(__ARM_NEON) && __ARM_NEON)
#		undef  MC_TARGET_HAVE_NEON
#		define MC_TARGET_HAVE_NEON 1
#		include <arm_neon.h>
#	endif

#	if defined(__APPLE__) && (!defined(__DARWIN_C_LEVEL) || __DARWIN_C_LEVEL < 199506L)
#		error "C99 compiler and Posix 1-2001 CRT required."
#	endif

#	undef MC_TARGET_APPLEXM
#	if    defined(__APPLE__)                              \
		&& defined(__IPHONE_OS_VERSION_MIN_ALLOWED)        \
		&& defined(__IPHONE_7_0)                           \
		&& __IPHONE_OS_VERSION_MIN_ALLOWED >= __IPHONE_7_0
#	define MC_TARGET_APPLEXM 1
#	endif

#	if    defined(__APPLE__)                                       \
		&& defined(MAC_OS_X_VERSION_MIN_REQUIRED)                   \
		&& defined(MAC_OS_X_VERSION_10_9)                           \
		&& MAC_OS_X_VERSION_MIN_REQUIRED   >= MAC_OS_X_VERSION_10_9
#	undef MC_TARGET_APPLEXM
#	define MC_TARGET_APPLEXM 1
#	endif

#	if defined(_MSC_VER)
#	include <intrin.h>
	MC_TARGET_FUNC int __builtin_ctz(unsigned int x)
	{
		unsigned long ret;
		unsigned __int32 m = mc_cast(unsigned __int32, x);
		_BitScanForward(&ret, m);
		return mc_cast(int, ret);
	}

	MC_TARGET_FUNC int __builtin_ctzll(unsigned long long x)
	{
		unsigned long ret;
		unsigned __int64 m = mc_cast(unsigned __int64, x);
		_BitScanForward64(&ret, m);
		return mc_cast(int, ret);
	}

	MC_TARGET_FUNC int __builtin_ctzl(unsigned long x)
	{
		return (sizeof(x) > sizeof(__int32) ?
			  __builtin_ctzll(mc_cast(unsigned long long, x))
			: __builtin_ctz(mc_cast(unsigned int, x))
		);
	}

	MC_TARGET_FUNC int __builtin_clz(unsigned int x)
	{
		unsigned long ret;
		unsigned __int32 m = mc_cast(unsigned __int32, x);
		_BitScanReverse(&ret, m);
		return mc_cast(int, ((CHAR_BIT * sizeof(unsigned __int32) - 1) ^ ret));
	}

	MC_TARGET_FUNC int __builtin_clzll(unsigned long long x)
	{
		unsigned long ret;
		unsigned __int64 m = mc_cast(unsigned __int64, x);
		_BitScanReverse64(&ret, x);
		return mc_cast(int, ((CHAR_BIT * sizeof(unsigned __int64) - 1) ^ ret));
	}

	MC_TARGET_FUNC int __builtin_clzl(unsigned long x)
	{
		unsigned long ret;
		if (sizeof(x) > sizeof(__int32)) {
			return __builtin_clzll(x);
		}
		unsigned __int32 m = mc_cast(unsigned __int32, x);
		_BitScanReverse(&ret, m);
		return mc_cast(int, ((CHAR_BIT * sizeof(unsigned __int32) - 1) ^ ret));
	}

	MC_TARGET_FUNC unsigned __int16 __builtin_bswap16(unsigned __int16 x)
	{ return _byteswap_ushort(x); }

	MC_TARGET_FUNC unsigned __int32 __builtin_bswap32(unsigned __int32 x)
	{ return _byteswap_ulong(x); }

	MC_TARGET_FUNC unsigned __int64 __builtin_bswap64(unsigned __int64 x)
	{ return _byteswap_uint64(x); }

#	undef __LITTLE_ENDIAN__
#	undef __BIG_ENDIAN__
#	undef __PDP_ENDIAN__
#	undef __ORDER_BIG_ENDIAN__
#	undef __ORDER_LITTLE_ENDIAN__
#	undef __ORDER_PDP_ENDIAN__
#	undef __BYTE_ORDER__

#	if defined(_M_PPC) || defined(_M_ALPHA)
#		define __BIG_ENDIAN__          1
#		define __ORDER_BIG_ENDIAN__    4321
#		define __ORDER_LITTLE_ENDIAN__ 1234
#		define __ORDER_PDP_ENDIAN__    3412
#		define __BYTE_ORDER__          __ORDER_BIG_ENDIAN__
#	elif  defined(__ATOM__)       \
		|| defined(__AVX__)        \
		|| defined(__AVX2__)       \
		|| defined(_M_IA64)        \
		|| defined(_M_I86)         \
		|| defined(_M_IX86)        \
		|| defined(_M_X64)         \
		|| defined(_M_AMD64)       \
		|| defined(_M_ARM)         \
		|| defined(_M_ARM_ARMV7VE) \
		|| defined(_M_ARM64)
#		define __LITTLE_ENDIAN__       1
#		define __ORDER_BIG_ENDIAN__    4321
#		define __ORDER_LITTLE_ENDIAN__ 1234
#		define __ORDER_PDP_ENDIAN__    3412
#		define __BYTE_ORDER__          __ORDER_LITTLE_ENDIAN__
#	else
#		error "Byteorder unknown."
#	endif
#	endif

#	if defined(__APPLE__) && defined(__MACH__)
#		include <machine/endian.h>
#	elif defined(__sun) || defined(sun) || defined(__SVR4)
#		include <sys/byteorder.h>
#	elif defined(__linux__)|| defined(__gnu_linux__)
#		include <endian.h>
#	elif defined(__unix__)      \
	|| defined(__bsdi__)        \
	|| defined(ANDROID)         \
	|| defined(__ANDROID__)     \
	|| defined(__ANDROID_API__) \
	|| defined(__FreeBSD__)     \
	|| defined(__NetBSD__)      \
	|| defined(__OpenBSD__)     \
	|| defined(__DragonFly__)
#		include <sys/endian.h>
#	else
#		error "Endian header not found."
#	endif

#	if !defined(BYTE_ORDER) && defined(__BYTE_ORDER__)
#		define BYTE_ORDER    __BYTE_ORDER__
#	endif
#	if !defined(BIG_ENDIAN) && defined(__ORDER_BIG_ENDIAN__)
#		define BIG_ENDIAN    __ORDER_BIG_ENDIAN__
#	endif
#	if !defined(LITTLE_ENDIAN) && defined(__ORDER_LITTLE_ENDIAN__)
#		define LITTLE_ENDIAN __ORDER_LITTLE_ENDIAN__
#	endif
#	if !defined(PDP_ENDIAN) && defined(__ORDER_PDP_ENDIAN__)
#		define PDP_ENDIAN    __ORDER_PDP_ENDIAN__
#	endif

#	if !defined(BYTE_ORDER)
#		error "Byteorder macro is not defined."
#	endif

#	define MC_TARGET_CTZ(x)     __builtin_ctz(x)
#	define MC_TARGET_CTZL(x)    __builtin_ctzl(x)
#	define MC_TARGET_CTZLL(x)   __builtin_ctzll(x)
#	define MC_TARGET_CLZ(x)     __builtin_clz(x)
#	define MC_TARGET_CLZL(x)    __builtin_clzl(x)
#	define MC_TARGET_CLZLL(x)   __builtin_clzll(x)

#	if (defined(__GNUC__) || defined(__clang__))
#		if ((__GNUC__ + 0 >= 4) || (__clang_major__ >= 4))
#			undef  MC_TARGET_SIGNBITF
#			undef  MC_TARGET_SIGNBIT
#			undef  MC_TARGET_SIGNBITL
#			undef  MC_TARGET_HAVE_SIGNBIT
#			define MC_TARGET_SIGNBITF(x) __builtin_signbitf(x)
#			define MC_TARGET_SIGNBIT(x)  __builtin_signbit(x)
#			define MC_TARGET_SIGNBITL(x) __builtin_signbitl(x)
#			define MC_TARGET_HAVE_SIGNBIT 1
#		else
#			if MC_TARGET_CCP98
#				if MC_TARGET_CCP11
#					define MC_TARGET_SIGNBITF(x) ::std::signbit(x)
#					define MC_TARGET_SIGNBIT(x)  ::std::signbit(x)
#					define MC_TARGET_SIGNBITL(x) ::std::signbit(x)
#				else
#					define MC_TARGET_SIGNBITF(x) ::signbit(x)
#					define MC_TARGET_SIGNBIT(x)  ::signbit(x)
#					define MC_TARGET_SIGNBITL(x) ::signbit(x)
#				endif
#			else
#				define MC_TARGET_SIGNBITF(x) signbit(x)
#				define MC_TARGET_SIGNBIT(x)  signbit(x)
#				define MC_TARGET_SIGNBITL(x) signbit(x)
#			endif
#			define MC_TARGET_HAVE_SIGNBIT 1
#		endif
#	else
#		if MC_TARGET_CCP98
#			if MC_TARGET_CCP11
#				define MC_TARGET_SIGNBITF(x) ::std::signbit(x)
#				define MC_TARGET_SIGNBIT(x)  ::std::signbit(x)
#				define MC_TARGET_SIGNBITL(x) ::std::signbit(x)
#			else
#				define MC_TARGET_SIGNBITF(x) ::signbit(x)
#				define MC_TARGET_SIGNBIT(x)  ::signbit(x)
#				define MC_TARGET_SIGNBITL(x) ::signbit(x)
#			endif
#		else
#			define MC_TARGET_SIGNBITF(x) signbit(x)
#			define MC_TARGET_SIGNBIT(x)  signbit(x)
#			define MC_TARGET_SIGNBITL(x) signbit(x)
#		endif
#		define MC_TARGET_HAVE_SIGNBIT 1
#	endif

#	define MC_TARGET_BSWAP16(x) __builtin_bswap16(x)
#	define MC_TARGET_BSWAP32(x) __builtin_bswap32(x)
#	define MC_TARGET_BSWAP64(x) __builtin_bswap64(x)

MC_TARGET_FUNC char               MC_TARGET_CHAR      (const char x)               { return x; }
MC_TARGET_FUNC short              MC_TARGET_SHORT     (const short x)              { return x; }
MC_TARGET_FUNC int                MC_TARGET_INT       (const int x)                { return x; }
MC_TARGET_FUNC long               MC_TARGET_LONG      (const long x)               { return x; }
MC_TARGET_FUNC unsigned char      MC_TARGET_UCHAR     (const unsigned char x)      { return x; }
MC_TARGET_FUNC unsigned short     MC_TARGET_USHORT    (const unsigned short x)     { return x; }
MC_TARGET_FUNC unsigned int       MC_TARGET_UINT      (const unsigned int x)       { return x; }
MC_TARGET_FUNC unsigned long      MC_TARGET_ULONG     (const unsigned long x)      { return x; }

#	if MC_TARGET_C99 || MC_TARGET_CPP11
MC_TARGET_FUNC long long          MC_TARGET_LONGLONG  (const long long x)          { return x; }
MC_TARGET_FUNC unsigned long long MC_TARGET_ULONGLONG (const unsigned long long x) { return x; }
#	else
#	define MC_TARGET_LONGLONG  MC_TARGET_LONG
#	define MC_TARGET_ULONGLONG MC_TARGET_ULONG
#	endif

#	ifndef MC_TARGET_ALLOCATOR_MAXSIZE
#		define MC_TARGET_ALLOCATOR_MAXSIZE UINT_MAX
#	endif

#endif /* !MC_TARGET_H */

/* EOF */
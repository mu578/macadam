//
// # -*- coding: utf-8, tab-width: 3 -*-

// mclimits.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_math.h>

#ifndef MCLIMITS_H
#define MCLIMITS_H

#pragma mark - numeric limits -

#	define MCLIMITS_EPSILONF MCK_KF(FLT_EPSILON)
#	define MCLIMITS_EPSILON  MCK_K(DBL_EPSILON)
#	define MCLIMITS_EPSILONL MCK_KL(LDBL_EPSILON)

#	define MCLIMITS_EXP_ARGMAXF    mc_cast_expr(const float   , mc_cast_expr(const int, FLT_MAX_EXP * M_LN2))
#	define MCLIMITS_EXP_ARGMAX     mc_cast_expr(const double  , mc_cast_expr(const int, DBL_MAX_EXP * M_LN2))
#	if (MC_TARGET_C99 || MC_TARGET_CPP17) && defined(M_LN2l)
#	define MCLIMITS_EXP_ARGMAXL mc_cast_expr(const long double, mc_cast_expr(const int, LDBL_MAX_EXP * M_LN2l))
#	else
#	define MCLIMITS_EXP_ARGMAXL mc_cast_expr(const long double, mc_cast_expr(const int, LDBL_MAX_EXP * M_LN2))
#	endif

#	define MCLIMITS_EXP_ARGMINF    mc_cast_expr(const float   , mc_cast_expr(const int, FLT_MIN_EXP * M_LN2))
#	define MCLIMITS_EXP_ARGMIN     mc_cast_expr(const double  , mc_cast_expr(const int, DBL_MIN_EXP * M_LN2))
#	if (MC_TARGET_C99 || MC_TARGET_CPP17) && defined(M_LN2l)
#	define MCLIMITS_EXP_ARGMINL mc_cast_expr(const long double, mc_cast_expr(const int, LDBL_MIN_EXP * M_LN2l))
#	else
#	define MCLIMITS_EXP_ARGMINL mc_cast_expr(const long double, mc_cast_expr(const int, LDBL_MIN_EXP * M_LN2))
#	endif

#	if MC_TARGET_CPP98
#	if MC_TARGET_CPP11
		static const float       MCLIMITS_LOWF  = ::std::sqrt(MCLIMITS_EPSILONF);
		static const double      MCLIMITS_LOW   = ::std::sqrt(MCLIMITS_EPSILON);
		static const long double MCLIMITS_LOWL  = ::std::sqrt(MCLIMITS_EPSILONL);

		static const float       MCLIMITS_TINYF = 1E-06f * MCLIMITS_LOWF;
		static const double      MCLIMITS_TINY  = 1E-06  * MCLIMITS_LOW;
		static const long double MCLIMITS_TINYL = 1E-06L * MCLIMITS_LOWL;
#	else
#		define MCLIMITS_LOWF  (mc_cast(const float      , ::sqrtf(MCLIMITS_EPSILONF)))
#		define MCLIMITS_LOW   (mc_cast(const double     , ::sqrt (MCLIMITS_EPSILON)))
#		define MCLIMITS_LOWL  (mc_cast(const long double, ::sqrtl(MCLIMITS_EPSILONL)))

#		define MCLIMITS_TINYF (mc_cast_expr(const float      , 1E-06f * MCLIMITS_LOWF))
#		define MCLIMITS_TINY  (mc_cast_expr(const double     , 1E-06  * MCLIMITS_LOW))
#		define MCLIMITS_TINYL (mc_cast_expr(const long double, 1E-06L * MCLIMITS_LOWL))
#	endif
#	else
#		define MCLIMITS_LOWF  (mc_cast(const float      , sqrtf(MCLIMITS_EPSILONF)))
#		define MCLIMITS_LOW   (mc_cast(const double     , sqrt (MCLIMITS_EPSILON)))
#		define MCLIMITS_LOWL  (mc_cast(const long double, sqrtl(MCLIMITS_EPSILONL)))

#		define MCLIMITS_TINYF (mc_cast_expr(const float      , 1E-06f * MCLIMITS_LOWF))
#		define MCLIMITS_TINY  (mc_cast_expr(const double     , 1E-06  * MCLIMITS_LOW))
#		define MCLIMITS_TINYL (mc_cast_expr(const long double, 1E-06L * MCLIMITS_LOWL))
#	endif

#	define MCLIMITS_MAXF     FLT_MAX
#	define MCLIMITS_MAX      DBL_MAX
#	define MCLIMITS_MAXL     LDBL_MAX

#	define MCLIMITS_MINF     FLT_MIN
#	define MCLIMITS_MIN      DBL_MIN
#	define MCLIMITS_MINL     LDBL_MIN

#	define MCLIMITS_BMAX     SCHAR_MAX
#	define MCLIMITS_SMAX     SHRT_MAX
#	define MCLIMITS_IMAX     INT_MAX
#	define MCLIMITS_LMAX     LONG_MAX
#	define MCLIMITS_LLMAX    LLONG_MAX

#	define MCLIMITS_UBMAX    UCHAR_MAX
#	define MCLIMITS_USMAX    USHRT_MAX
#	define MCLIMITS_UIMAX    UINT_MAX
#	define MCLIMITS_ULMAX    ULONG_MAX
#	define MCLIMITS_ULLMAX   ULLONG_MAX

#	define MCLIMITS_BMIN     SCHAR_MIN
#	define MCLIMITS_SMIN     SHRT_MIN
#	define MCLIMITS_IMIN     INT_MIN
#	define MCLIMITS_LMIN     LONG_MIN
#	define MCLIMITS_LLMIN    LLONG_MIN

#	define MCLIMITS_UBMIN    (0)
#	define MCLIMITS_USMIN    (0)
#	define MCLIMITS_UIMIN    (0)
#	define MCLIMITS_ULMIN    (0)
#	define MCLIMITS_ULLMIN   (0)

#	define MCLIMITS_CBITS    CHAR_BIT

#pragma mark - mclimits_epsilonof -

#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T           mclimits_epsilonof              (const T& x)           { mc_unused(x); return 0;                 }
template <>        MC_TARGET_INLINE float       mclimits_epsilonof<float>       (const float& x)       { mc_unused(x); return MCLIMITS_EPSILONF; }
template <>        MC_TARGET_INLINE double      mclimits_epsilonof<double>      (const double& x)      { mc_unused(x); return MCLIMITS_EPSILON;  }
template <>        MC_TARGET_INLINE long double mclimits_epsilonof<long double> (const long double& x) { mc_unused(x); return MCLIMITS_EPSILONL; }
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float       mclimits_epsilonof (const float x)       { mc_unused(x); return MCLIMITS_EPSILONF; }
MC_TARGET_ALIAS double      mclimits_epsilonof (const double x)      { mc_unused(x); return MCLIMITS_EPSILON;  }
MC_TARGET_ALIAS long double mclimits_epsilonof (const long double x) { mc_unused(x); return MCLIMITS_EPSILONL; }
#	elif MC_TARGET_C11
#	define mclimits_epsilonof(x) _Generic(x  \
	, float              : MCLIMITS_EPSILONF \
	, double             : MCLIMITS_EPSILON  \
	, long double        : MCLIMITS_EPSILONL \
	, signed char        : (0)               \
	, short              : (0)               \
	, int                : (0)               \
	, long               : (0)               \
	, long long          : (0)               \
	, unsigned char      : (0)               \
	, unsigned short     : (0)               \
	, unsigned int       : (0)               \
	, unsigned long      : (0)               \
	, unsigned long long : (0)               \
)
#	elif MC_TARGET_HAVE_TYPEOF
#	define mclimits_epsilonof(x) mc_cast(MC_TARGET_TYPEOF(x),                     \
	(                                                                             \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)       ? MCLIMITS_EPSILONF \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)      ? MCLIMITS_EPSILON  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double) ? MCLIMITS_EPSILONL \
		: 0                                                                        \
	))
#	else
#	define mclimits_epsilonof(x) (0)
#	endif

#pragma mark - mclimits_maxof -

#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T                  mclimits_maxof                     (const T& x)                  { mc_unused(x); return 0;               }
template <>        MC_TARGET_INLINE float              mclimits_maxof<float>              (const float& x)              { mc_unused(x); return MCLIMITS_MAXF;   }
template <>        MC_TARGET_INLINE double             mclimits_maxof<double>             (const double& x)             { mc_unused(x); return MCLIMITS_MAX;    }
template <>        MC_TARGET_INLINE long double        mclimits_maxof<long double>        (const long double& x)        { mc_unused(x); return MCLIMITS_MAXL;   }
template <>        MC_TARGET_INLINE signed char        mclimits_maxof<signed char>        (const signed char& x)        { mc_unused(x); return MCLIMITS_BMAX;   }
template <>        MC_TARGET_INLINE short              mclimits_maxof<short>              (const short& x)              { mc_unused(x); return MCLIMITS_SMAX;   }
template <>        MC_TARGET_INLINE int                mclimits_maxof<int>                (const int& x)                { mc_unused(x); return MCLIMITS_IMAX;   }
template <>        MC_TARGET_INLINE long               mclimits_maxof<long>               (const long& x)               { mc_unused(x); return MCLIMITS_LMAX;   }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE long long          mclimits_maxof<long long>          (const long long& x)          { mc_unused(x); return MCLIMITS_LLMAX;  }
#	endif
template <>        MC_TARGET_INLINE unsigned char      mclimits_maxof<unsigned char>      (const unsigned char& x)      { mc_unused(x); return MCLIMITS_UBMAX;  }
template <>        MC_TARGET_INLINE unsigned short     mclimits_maxof<unsigned short>     (const unsigned short& x)     { mc_unused(x); return MCLIMITS_USMAX;  }
template <>        MC_TARGET_INLINE unsigned int       mclimits_maxof<unsigned int>       (const unsigned int& x)       { mc_unused(x); return MCLIMITS_UIMAX;  }
template <>        MC_TARGET_INLINE unsigned long      mclimits_maxof<unsigned long>      (const unsigned long& x)      { mc_unused(x); return MCLIMITS_ULMAX;  }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE unsigned long long mclimits_maxof<unsigned long long> (const unsigned long long& x) { mc_unused(x); return MCLIMITS_ULLMAX; }
#	endif
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float              mclimits_maxof (const float x)              { mc_unused(x); return MCLIMITS_MAXF;   }
MC_TARGET_ALIAS double             mclimits_maxof (const double x)             { mc_unused(x); return MCLIMITS_MAX;    }
MC_TARGET_ALIAS long double        mclimits_maxof (const long double x)        { mc_unused(x); return MCLIMITS_MAXL;   }
MC_TARGET_ALIAS signed char        mclimits_maxof (const signed char x)        { mc_unused(x); return MCLIMITS_BMAX;   }
MC_TARGET_ALIAS short              mclimits_maxof (const short x)              { mc_unused(x); return MCLIMITS_SMAX;   }
MC_TARGET_ALIAS int                mclimits_maxof (const int x)                { mc_unused(x); return MCLIMITS_IMAX;   }
MC_TARGET_ALIAS long               mclimits_maxof (const long x)               { mc_unused(x); return MCLIMITS_LMAX;   }
#	if MC_TARGET_C99
MC_TARGET_ALIAS long long          mclimits_maxof (const long long x)          { mc_unused(x); return MCLIMITS_LLMAX;  }
#	endif
MC_TARGET_ALIAS unsigned char      mclimits_maxof (const unsigned char x)      { mc_unused(x); return MCLIMITS_UBMAX;  }
MC_TARGET_ALIAS unsigned short     mclimits_maxof (const unsigned short x)     { mc_unused(x); return MCLIMITS_USMAX;  }
MC_TARGET_ALIAS unsigned int       mclimits_maxof (const unsigned int x)       { mc_unused(x); return MCLIMITS_UIMAX;  }
MC_TARGET_ALIAS unsigned long      mclimits_maxof (const unsigned long x)      { mc_unused(x); return MCLIMITS_ULMAX;  }
#	if MC_TARGET_C99
MC_TARGET_ALIAS unsigned long long mclimits_maxof (const unsigned long long x) { mc_unused(x); return MCLIMITS_ULLMAX; }
#	endif
#	elif MC_TARGET_C11
#	define mclimits_maxof(x) _Generic(x    \
	, float              : MCLIMITS_MAXF   \
	, double             : MCLIMITS_MAX    \
	, long double        : MCLIMITS_MAXL   \
	, signed char        : MCLIMITS_BMAX   \
	, short              : MCLIMITS_SMAX   \
	, int                : MCLIMITS_IMAX   \
	, long               : MCLIMITS_LMAX   \
	, long long          : MCLIMITS_LLMAX  \
	, unsigned char      : MCLIMITS_UBMAX  \
	, unsigned short     : MCLIMITS_USMAX  \
	, unsigned int       : MCLIMITS_UIMAX  \
	, unsigned long      : MCLIMITS_ULMAX  \
	, unsigned long long : MCLIMITS_ULLMAX \
)
#	elif MC_TARGET_HAVE_TYPEOF
#	if MC_TARGET_C99
#	define mclimits_maxof(x) mc_cast(MC_TARGET_TYPEOF(x),                              \
	(                                                                                  \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? MCLIMITS_MAXF   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? MCLIMITS_MAX    \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? MCLIMITS_MAXL   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? MCLIMITS_BMAX   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? MCLIMITS_SMAX   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? MCLIMITS_IMAX   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? MCLIMITS_LLMAX  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long long)          ? MCLIMITS_LLMAX  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? MCLIMITS_UBMAX  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? MCLIMITS_USMAX  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? MCLIMITS_UIMAX  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? MCLIMITS_ULLMAX \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long long) ? MCLIMITS_ULLMAX \
		: 0                                                                             \
	))
#	else
#	define mclimits_maxof(x) mc_cast(MC_TARGET_TYPEOF(x),                              \
	(                                                                                  \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? MCLIMITS_MAXF   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? MCLIMITS_MAX    \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? MCLIMITS_MAXL   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? MCLIMITS_BMAX   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? MCLIMITS_SMAX   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? MCLIMITS_IMAX   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? MCLIMITS_LLMAX  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? MCLIMITS_UBMAX  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? MCLIMITS_USMAX  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? MCLIMITS_UIMAX  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? MCLIMITS_ULLMAX \
		: 0                                                                             \
	))
#	endif
#	else
#	define mclimits_maxof(x) (0)
#	endif

#pragma mark - mclimits_minof -

#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T                  mclimits_minof                     (const T& x)                  { mc_unused(x); return 0;               }
template <>        MC_TARGET_INLINE float              mclimits_minof<float>              (const float& x)              { mc_unused(x); return MCLIMITS_MINF;   }
template <>        MC_TARGET_INLINE double             mclimits_minof<double>             (const double& x)             { mc_unused(x); return MCLIMITS_MIN;    }
template <>        MC_TARGET_INLINE long double        mclimits_minof<long double>        (const long double& x)        { mc_unused(x); return MCLIMITS_MINL;   }
template <>        MC_TARGET_INLINE signed char        mclimits_minof<signed char>        (const signed char& x)        { mc_unused(x); return MCLIMITS_BMIN;   }
template <>        MC_TARGET_INLINE short              mclimits_minof<short>              (const short& x)              { mc_unused(x); return MCLIMITS_SMIN;   }
template <>        MC_TARGET_INLINE int                mclimits_minof<int>                (const int& x)                { mc_unused(x); return MCLIMITS_IMIN;   }
template <>        MC_TARGET_INLINE long               mclimits_minof<long>               (const long& x)               { mc_unused(x); return MCLIMITS_LMIN;   }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE long long          mclimits_minof<long long>          (const long long& x)          { mc_unused(x); return MCLIMITS_LLMIN;  }
#	endif

template <>        MC_TARGET_INLINE unsigned char      mclimits_minof<unsigned char>      (const unsigned char& x)      { mc_unused(x); return MCLIMITS_UBMIN;  }
template <>        MC_TARGET_INLINE unsigned short     mclimits_minof<unsigned short>     (const unsigned short& x)     { mc_unused(x); return MCLIMITS_USMIN;  }
template <>        MC_TARGET_INLINE unsigned int       mclimits_minof<unsigned int>       (const unsigned int& x)       { mc_unused(x); return MCLIMITS_UIMIN;  }
template <>        MC_TARGET_INLINE unsigned long      mclimits_minof<unsigned long>      (const unsigned long& x)      { mc_unused(x); return MCLIMITS_ULMIN;  }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE unsigned long long mclimits_minof<unsigned long long> (const unsigned long long& x) { mc_unused(x); return MCLIMITS_ULLMIN; }
#	endif
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float              mclimits_minof (const float x)              { mc_unused(x); return MCLIMITS_MINF;   }
MC_TARGET_ALIAS double             mclimits_minof (const double x)             { mc_unused(x); return MCLIMITS_MIN;    }
MC_TARGET_ALIAS long double        mclimits_minof (const long double x)        { mc_unused(x); return MCLIMITS_MINL;   }
MC_TARGET_ALIAS signed char        mclimits_minof (const signed char x)        { mc_unused(x); return MCLIMITS_BMIN;   }
MC_TARGET_ALIAS short              mclimits_minof (const short x)              { mc_unused(x); return MCLIMITS_SMIN;   }
MC_TARGET_ALIAS int                mclimits_minof (const int x)                { mc_unused(x); return MCLIMITS_IMIN;   }
MC_TARGET_ALIAS long               mclimits_minof (const long x)               { mc_unused(x); return MCLIMITS_LMIN;   }
#	if MC_TARGET_C99
MC_TARGET_ALIAS long long          mclimits_minof (const long long x)          { mc_unused(x); return MCLIMITS_LLMIN;  }
#	endif
MC_TARGET_ALIAS unsigned char      mclimits_minof (const unsigned char x)      { mc_unused(x); return MCLIMITS_UBMIN;  }
MC_TARGET_ALIAS unsigned short     mclimits_minof (const unsigned short x)     { mc_unused(x); return MCLIMITS_USMIN;  }
MC_TARGET_ALIAS unsigned int       mclimits_minof (const unsigned int x)       { mc_unused(x); return MCLIMITS_UIMIN;  }
MC_TARGET_ALIAS unsigned long      mclimits_minof (const unsigned long x)      { mc_unused(x); return MCLIMITS_ULMIN;  }
#	if MC_TARGET_C99
MC_TARGET_ALIAS unsigned long long mclimits_minof (const unsigned long long x) { mc_unused(x); return MCLIMITS_ULLMIN; }
#	endif
#	elif MC_TARGET_C11
#	define mclimits_minof(x) _Generic(x    \
	, float              : MCLIMITS_MINF   \
	, double             : MCLIMITS_MIN    \
	, long double        : MCLIMITS_MINL   \
	, signed char        : MCLIMITS_BMIN   \
	, short              : MCLIMITS_SMIN   \
	, int                : MCLIMITS_IMIN   \
	, long               : MCLIMITS_LMIN   \
	, long long          : MCLIMITS_LLMIN  \
	, unsigned char      : MCLIMITS_UBMIN  \
	, unsigned short     : MCLIMITS_USMIN  \
	, unsigned int       : MCLIMITS_UIMIN  \
	, unsigned long      : MCLIMITS_ULMIN  \
	, unsigned long long : MCLIMITS_ULLMIN \
)
#	elif MC_TARGET_HAVE_TYPEOF
#	if MC_TARGET_C99
#	define mclimits_minof(x) mc_cast(MC_TARGET_TYPEOF(x),                              \
	(                                                                                  \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? MCLIMITS_MINF   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? MCLIMITS_MIN    \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? MCLIMITS_MINL   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? MCLIMITS_BMIN   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? MCLIMITS_SMIN   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? MCLIMITS_IMIN   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? MCLIMITS_LLMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long long)          ? MCLIMITS_LLMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? MCLIMITS_UBMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? MCLIMITS_USMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? MCLIMITS_UIMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? MCLIMITS_ULLMIN \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long long) ? MCLIMITS_ULLMIN \
		: 0                                                                             \
	))
#	else
#	define mclimits_minof(x) mc_cast(MC_TARGET_TYPEOF(x),                              \
	(                                                                                  \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? MCLIMITS_MINF   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? MCLIMITS_MIN    \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? MCLIMITS_MINL   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? MCLIMITS_BMIN   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? MCLIMITS_SMIN   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? MCLIMITS_IMIN   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? MCLIMITS_LLMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? MCLIMITS_UBMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? MCLIMITS_USMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? MCLIMITS_UIMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? MCLIMITS_ULLMIN \
		: 0                                                                             \
	))
#	endif
#	else
#	define mclimits_minof(x) (0)
#	endif

#pragma mark - mclimits_lowestof -

#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE T                  mclimits_lowestof                     (const T& x)                  { mc_unused(x); return  0;               }
template <>        MC_TARGET_INLINE float              mclimits_lowestof<float>              (const float& x)              { mc_unused(x); return -MCLIMITS_MAXF;   }
template <>        MC_TARGET_INLINE double             mclimits_lowestof<double>             (const double& x)             { mc_unused(x); return -MCLIMITS_MAX;    }
template <>        MC_TARGET_INLINE long double        mclimits_lowestof<long double>        (const long double& x)        { mc_unused(x); return -MCLIMITS_MAXL;   }
template <>        MC_TARGET_INLINE signed char        mclimits_lowestof<signed char>        (const signed char& x)        { mc_unused(x); return  MCLIMITS_BMIN;   }
template <>        MC_TARGET_INLINE short              mclimits_lowestof<short>              (const short& x)              { mc_unused(x); return  MCLIMITS_SMIN;   }
template <>        MC_TARGET_INLINE int                mclimits_lowestof<int>                (const int& x)                { mc_unused(x); return  MCLIMITS_IMIN;   }
template <>        MC_TARGET_INLINE long               mclimits_lowestof<long>               (const long& x)               { mc_unused(x); return  MCLIMITS_LMIN;   }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE long long          mclimits_lowestof<long long>          (const long long& x)          { mc_unused(x); return  MCLIMITS_LLMIN;  }
#	endif
template <>        MC_TARGET_INLINE unsigned char      mclimits_lowestof<unsigned char>      (const unsigned char& x)      { mc_unused(x); return  MCLIMITS_UBMIN;  }
template <>        MC_TARGET_INLINE unsigned short     mclimits_lowestof<unsigned short>     (const unsigned short& x)     { mc_unused(x); return  MCLIMITS_USMIN;  }
template <>        MC_TARGET_INLINE unsigned int       mclimits_lowestof<unsigned int>       (const unsigned int& x)       { mc_unused(x); return  MCLIMITS_UIMIN;  }
template <>        MC_TARGET_INLINE unsigned long      mclimits_lowestof<unsigned long>      (const unsigned long& x)      { mc_unused(x); return  MCLIMITS_ULMIN;  }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE unsigned long long mclimits_lowestof<unsigned long long> (const unsigned long long& x) { mc_unused(x); return  MCLIMITS_ULLMIN; }
#	endif
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS float              mclimits_lowestof (const float x)              { mc_unused(x); return -MCLIMITS_MAXF;   }
MC_TARGET_ALIAS double             mclimits_lowestof (const double x)             { mc_unused(x); return -MCLIMITS_MAX;    }
MC_TARGET_ALIAS long double        mclimits_lowestof (const long double x)        { mc_unused(x); return -MCLIMITS_MAXL;   }
MC_TARGET_ALIAS signed char        mclimits_lowestof (const signed char x)        { mc_unused(x); return  MCLIMITS_BMIN;   }
MC_TARGET_ALIAS short              mclimits_lowestof (const short x)              { mc_unused(x); return  MCLIMITS_SMIN;   }
MC_TARGET_ALIAS int                mclimits_lowestof (const int x)                { mc_unused(x); return  MCLIMITS_IMIN;   }
MC_TARGET_ALIAS long               mclimits_lowestof (const long x)               { mc_unused(x); return  MCLIMITS_LMIN;   }
#	if MC_TARGET_C99
MC_TARGET_ALIAS long long          mclimits_lowestof (const long long x)          { mc_unused(x); return  MCLIMITS_LLMIN;  }
#	endif
MC_TARGET_ALIAS unsigned char      mclimits_lowestof (const unsigned char x)      { mc_unused(x); return  MCLIMITS_UBMIN;  }
MC_TARGET_ALIAS unsigned short     mclimits_lowestof (const unsigned short x)     { mc_unused(x); return  MCLIMITS_USMIN;  }
MC_TARGET_ALIAS unsigned int       mclimits_lowestof (const unsigned int x)       { mc_unused(x); return  MCLIMITS_UIMIN;  }
MC_TARGET_ALIAS unsigned long      mclimits_lowestof (const unsigned long x)      { mc_unused(x); return  MCLIMITS_ULMIN;  }
#	if MC_TARGET_C99
MC_TARGET_ALIAS unsigned long long mclimits_lowestof (const unsigned long long x) { mc_unused(x); return  MCLIMITS_ULLMIN; }
#	endif
#	elif MC_TARGET_C11
#	define mclimits_lowestof(x) _Generic(x  \
	, float              : -MCLIMITS_MAXF   \
	, double             : -MCLIMITS_MAX    \
	, long double        : -MCLIMITS_MAXL   \
	, signed char        :  MCLIMITS_BMIN   \
	, short              :  MCLIMITS_SMIN   \
	, int                :  MCLIMITS_IMIN   \
	, long               :  MCLIMITS_LMIN   \
	, long long          :  MCLIMITS_LLMIN  \
	, unsigned char      :  MCLIMITS_UBMIN  \
	, unsigned short     :  MCLIMITS_USMIN  \
	, unsigned int       :  MCLIMITS_UIMIN  \
	, unsigned long      :  MCLIMITS_ULMIN  \
	, unsigned long long :  MCLIMITS_ULLMIN \
)
#	elif MC_TARGET_HAVE_TYPEOF
#	if MC_TARGET_C99
#	define mclimits_lowestof(x) mc_cast(MC_TARGET_TYPEOF(x),                            \
	(                                                                                   \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? -MCLIMITS_MAXF   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? -MCLIMITS_MAX    \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? -MCLIMITS_MAXL   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ?  MCLIMITS_BMIN   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ?  MCLIMITS_SMIN   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ?  MCLIMITS_IMIN   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ?  MCLIMITS_LLMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long long)          ?  MCLIMITS_LLMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ?  MCLIMITS_UBMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ?  MCLIMITS_USMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ?  MCLIMITS_UIMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ?  MCLIMITS_ULLMIN \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long long) ?  MCLIMITS_ULLMIN \
		: 0                                                                              \
	))
#	else
#	define mclimits_lowestof(x) mc_cast(MC_TARGET_TYPEOF(x),                            \
	(                                                                                   \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? -MCLIMITS_MAXF   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? -MCLIMITS_MAX    \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? -MCLIMITS_MAXL   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ?  MCLIMITS_BMIN   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ?  MCLIMITS_SMIN   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ?  MCLIMITS_IMIN   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ?  MCLIMITS_LLMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ?  MCLIMITS_UBMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ?  MCLIMITS_USMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ?  MCLIMITS_UIMIN  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ?  MCLIMITS_ULLMIN \
		: 0                                                                              \
	))
#	endif
#	else
#	define mclimits_lowestof(x) (0)
#	endif

#pragma mark - mclimits_maxexpof -

#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE int mclimits_maxexpof                     (const T& x)                  { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_maxexpof<float>              (const float& x)              { mc_unused(x); return mc_cast(int, FLT_MAX_EXP);  }
template <>        MC_TARGET_INLINE int mclimits_maxexpof<double>             (const double& x)             { mc_unused(x); return mc_cast(int, DBL_MAX_EXP);  }
template <>        MC_TARGET_INLINE int mclimits_maxexpof<long double>        (const long double& x)        { mc_unused(x); return mc_cast(int, LDBL_MAX_EXP); }
template <>        MC_TARGET_INLINE int mclimits_maxexpof<signed char>        (const signed char& x)        { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_maxexpof<short>              (const short& x)              { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_maxexpof<int>                (const int& x)                { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_maxexpof<long>               (const long& x)               { mc_unused(x); return 0;                          }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE int mclimits_maxexpof<long long>          (const long long& x)          { mc_unused(x); return 0;                          }
#	endif
template <>        MC_TARGET_INLINE int mclimits_maxexpof<unsigned char>      (const unsigned char& x)      { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_maxexpof<unsigned short>     (const unsigned short& x)     { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_maxexpof<unsigned int>       (const unsigned int& x)       { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_maxexpof<unsigned long>      (const unsigned long& x)      { mc_unused(x); return 0;                          }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE int mclimits_maxexpof<unsigned long long> (const unsigned long long& x) { mc_unused(x); return 0;                          }
#	endif
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS int mclimits_maxexpof (const float x)              { mc_unused(x); return mc_cast(int, FLT_MAX_EXP);  }
MC_TARGET_ALIAS int mclimits_maxexpof (const double x)             { mc_unused(x); return mc_cast(int, DBL_MAX_EXP);  }
MC_TARGET_ALIAS int mclimits_maxexpof (const long double x)        { mc_unused(x); return mc_cast(int, LDBL_MAX_EXP); }
MC_TARGET_ALIAS int mclimits_maxexpof (const signed char x)        { mc_unused(x); return 0;                          }
MC_TARGET_ALIAS int mclimits_maxexpof (const short x)              { mc_unused(x); return 0;                          }
MC_TARGET_ALIAS int mclimits_maxexpof (const int x)                { mc_unused(x); return 0;                          }
MC_TARGET_ALIAS int mclimits_maxexpof (const long x)               { mc_unused(x); return 0;                          }
#	if MC_TARGET_C99
MC_TARGET_ALIAS int mclimits_maxexpof (const long long x)          { mc_unused(x); return 0;                          }
#	endif
MC_TARGET_ALIAS int mclimits_maxexpof (const unsigned char x)      { mc_unused(x); return 0;                          }
MC_TARGET_ALIAS int mclimits_maxexpof (const unsigned short x)     { mc_unused(x); return 0;                          }
MC_TARGET_ALIAS int mclimits_maxexpof (const unsigned int x)       { mc_unused(x); return 0;                          }
MC_TARGET_ALIAS int mclimits_maxexpof (const unsigned long x)      { mc_unused(x); return 0;                          }
#	if MC_TARGET_C99
MC_TARGET_ALIAS int mclimits_maxexpof (const unsigned long long x) { mc_unused(x); return 0;                          }
#	endif
#	elif MC_TARGET_C11
#	define mclimits_maxexpof(x) mc_cast(int, _Generic(x \
	, float              : FLT_MAX_EXP                  \
	, double             : DBL_MAX_EXP                  \
	, long double        : LDBL_MAX_EXP                 \
	, signed char        : 0                            \
	, short              : 0                            \
	, int                : 0                            \
	, long               : 0                            \
	, long long          : 0                            \
	, unsigned char      : 0                            \
	, unsigned short     : 0                            \
	, unsigned int       : 0                            \
	, unsigned long      : 0                            \
	, unsigned long long : 0                            \
))
#	elif MC_TARGET_HAVE_TYPEOF
#	if MC_TARGET_C99
#	define mclimits_maxexpof(x) mc_cast(int,                                        \
	(                                                                               \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? FLT_MAX_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? DBL_MAX_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? LDBL_MAX_EXP \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long long)          ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long long) ? 0            \
		: 0                                                                          \
	))
#	else
#	define mclimits_maxexpof(x) mc_cast(int,                                        \
	(                                                                               \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? FLT_MAX_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? DBL_MAX_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? LDBL_MAX_EXP \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? 0            \
		: 0                                                                          \
	))
#	endif
#	else
#	define mclimits_maxexpof(x) (0)
#	endif

#pragma mark - mclimits_minexpof -
#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE int mclimits_minexpof                     (const T& x)                  { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_minexpof<float>              (const float& x)              { mc_unused(x); return mc_cast(int, FLT_MIN_EXP);  }
template <>        MC_TARGET_INLINE int mclimits_minexpof<double>             (const double& x)             { mc_unused(x); return mc_cast(int, DBL_MIN_EXP);  }
template <>        MC_TARGET_INLINE int mclimits_minexpof<long double>        (const long double& x)        { mc_unused(x); return mc_cast(int, LDBL_MIN_EXP); }
template <>        MC_TARGET_INLINE int mclimits_minexpof<signed char>        (const signed char& x)        { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_minexpof<short>              (const short& x)              { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_minexpof<int>                (const int& x)                { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_minexpof<long>               (const long& x)               { mc_unused(x); return 0;                          }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE int mclimits_minexpof<long long>          (const long long& x)          { mc_unused(x); return 0;                          }
#	endif
template <>        MC_TARGET_INLINE int mclimits_minexpof<unsigned char>      (const unsigned char& x)      { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_minexpof<unsigned short>     (const unsigned short& x)     { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_minexpof<unsigned int>       (const unsigned int& x)       { mc_unused(x); return 0;                          }
template <>        MC_TARGET_INLINE int mclimits_minexpof<unsigned long>      (const unsigned long& x)      { mc_unused(x); return 0;                          }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE int mclimits_minexpof<unsigned long long> (const unsigned long long& x) { mc_unused(x); return 0;                          }
#	endif
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS int mclimits_minexpof (const float x)              { mc_unused(x); return mc_cast(int, FLT_MIN_EXP);  }
MC_TARGET_ALIAS int mclimits_minexpof (const double x)             { mc_unused(x); return mc_cast(int, DBL_MIN_EXP);  }
MC_TARGET_ALIAS int mclimits_minexpof (const long double x)        { mc_unused(x); return mc_cast(int, LDBL_MIN_EXP); }
MC_TARGET_ALIAS int mclimits_minexpof (const signed char x)        { mc_unused(x); return 0;                          }
MC_TARGET_ALIAS int mclimits_minexpof (const short x)              { mc_unused(x); return 0;                          }
MC_TARGET_ALIAS int mclimits_minexpof (const int x)                { mc_unused(x); return 0;                          }
MC_TARGET_ALIAS int mclimits_minexpof (const long x)               { mc_unused(x); return 0;                          }
#	if MC_TARGET_C99
MC_TARGET_ALIAS int mclimits_minexpof (const long long x)          { mc_unused(x); return 0;                          }
#	endif
MC_TARGET_ALIAS int mclimits_minexpof (const unsigned char x)      { mc_unused(x); return 0;                          }
MC_TARGET_ALIAS int mclimits_minexpof (const unsigned short x)     { mc_unused(x); return 0;                          }
MC_TARGET_ALIAS int mclimits_minexpof (const unsigned int x)       { mc_unused(x); return 0;                          }
MC_TARGET_ALIAS int mclimits_minexpof (const unsigned long x)      { mc_unused(x); return 0;                          }
#	if MC_TARGET_C99
MC_TARGET_ALIAS int mclimits_minexpof (const unsigned long long x) { mc_unused(x); return 0;                          }
#	endif
#	elif MC_TARGET_C11
#	define mclimits_minexpof(x) mc_cast(int, _Generic(x \
	, float              : FLT_MIN_EXP                  \
	, double             : DBL_MIN_EXP                  \
	, long double        : LDBL_MIN_EXP                 \
	, signed char        : 0                            \
	, short              : 0                            \
	, int                : 0                            \
	, long               : 0                            \
	, long long          : 0                            \
	, unsigned char      : 0                            \
	, unsigned short     : 0                            \
	, unsigned int       : 0                            \
	, unsigned long      : 0                            \
	, unsigned long long : 0                            \
))s
#	elif MC_TARGET_HAVE_TYPEOF
#	if MC_TARGET_C99
#	define mclimits_minexpof(x) mc_cast(int,                                        \
	(                                                                               \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? FLT_MIN_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? DBL_MIN_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? LDBL_MIN_EXP \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long long)          ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long long) ? 0            \
		: 0                                                                          \
	))
#	else
#	define mclimits_minexpof(x) mc_cast(int,                                        \
	(                                                                               \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? FLT_MIN_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? DBL_MIN_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? LDBL_MIN_EXP \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? 0            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? 0            \
		: 0                                                                          \
	))
#	endif
#	else
#	define mclimits_minexpof(x) (0)
#	endif

#pragma mark - mclimits_maxexp10of -

#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE int mclimits_maxexp10of                     (const T& x)                  { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<float>              (const float& x)              { mc_unused(x); return mc_cast(int, FLT_MAX_10_EXP);  }
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<double>             (const double& x)             { mc_unused(x); return mc_cast(int, DBL_MAX_10_EXP);  }
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<long double>        (const long double& x)        { mc_unused(x); return mc_cast(int, LDBL_MAX_10_EXP); }
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<signed char>        (const signed char& x)        { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<short>              (const short& x)              { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<int>                (const int& x)                { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<long>               (const long& x)               { mc_unused(x); return 0;                             }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<long long>          (const long long& x)          { mc_unused(x); return 0;                             }
#	endif
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<unsigned char>      (const unsigned char& x)      { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<unsigned short>     (const unsigned short& x)     { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<unsigned int>       (const unsigned int& x)       { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<unsigned long>      (const unsigned long& x)      { mc_unused(x); return 0;                             }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE int mclimits_maxexp10of<unsigned long long> (const unsigned long long& x) { mc_unused(x); return 0;                             }
#	endif
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS int mclimits_maxexp10of (const float x)              { mc_unused(x); return mc_cast(int, FLT_MAX_10_EXP);  }
MC_TARGET_ALIAS int mclimits_maxexp10of (const double x)             { mc_unused(x); return mc_cast(int, DBL_MAX_10_EXP);  }
MC_TARGET_ALIAS int mclimits_maxexp10of (const long double x)        { mc_unused(x); return mc_cast(int, LDBL_MAX_10_EXP); }
MC_TARGET_ALIAS int mclimits_maxexp10of (const signed char x)        { mc_unused(x); return 0;                             }
MC_TARGET_ALIAS int mclimits_maxexp10of (const short x)              { mc_unused(x); return 0;                             }
MC_TARGET_ALIAS int mclimits_maxexp10of (const int x)                { mc_unused(x); return 0;                             }
MC_TARGET_ALIAS int mclimits_maxexp10of (const long x)               { mc_unused(x); return 0;                             }
#	if MC_TARGET_C99
MC_TARGET_ALIAS int mclimits_maxexp10of (const long long x)          { mc_unused(x); return 0;                             }
#	endif
MC_TARGET_ALIAS int mclimits_maxexp10of (const unsigned char x)      { mc_unused(x); return 0;                             }
MC_TARGET_ALIAS int mclimits_maxexp10of (const unsigned short x)     { mc_unused(x); return 0;                             }
MC_TARGET_ALIAS int mclimits_maxexp10of (const unsigned int x)       { mc_unused(x); return 0;                             }
MC_TARGET_ALIAS int mclimits_maxexp10of (const unsigned long x)      { mc_unused(x); return 0;                             }
#	if MC_TARGET_C99
MC_TARGET_ALIAS int mclimits_maxexp10of (const unsigned long long x) { mc_unused(x); return 0;                             }
#	endif
#	elif MC_TARGET_C11
#	define mclimits_maxexp10of(x) mc_cast(int, _Generic(x \
	, float              : FLT_MAX_10_EXP                 \
	, double             : DBL_MAX_10_EXP                 \
	, long double        : LDBL_MAX_10_EXP                \
	, signed char        : 0                              \
	, short              : 0                              \
	, int                : 0                              \
	, long               : 0                              \
	, long long          : 0                              \
	, unsigned char      : 0                              \
	, unsigned short     : 0                              \
	, unsigned int       : 0                              \
	, unsigned long      : 0                              \
	, unsigned long long : 0                              \
))
#	elif MC_TARGET_HAVE_TYPEOF
#	if MC_TARGET_C99
#	define mclimits_maxexp10of(x) mc_cast(int,                                         \
	(                                                                                  \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? FLT_MAX_10_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? DBL_MAX_10_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? LDBL_MAX_10_EXP \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long long)          ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long long) ? 0               \
		: 0                                                                             \
	))
#	else
#	define mclimits_maxexp10of(x) mc_cast(int,                                         \
	(                                                                                  \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? FLT_MAX_10_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? DBL_MAX_10_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? LDBL_MAX_10_EXP \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? 0               \
		: 0                                                                             \
	))
#	endif
#	else
#	define mclimits_maxexp10of(x) (0)
#	endif

#pragma mark - mclimits_minexp10of -

#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE int mclimits_minexp10of                     (const T& x)                  { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_minexp10of<float>              (const float& x)              { mc_unused(x); return mc_cast(int, FLT_MIN_10_EXP);  }
template <>        MC_TARGET_INLINE int mclimits_minexp10of<double>             (const double& x)             { mc_unused(x); return mc_cast(int, DBL_MIN_10_EXP);  }
template <>        MC_TARGET_INLINE int mclimits_minexp10of<long double>        (const long double& x)        { mc_unused(x); return mc_cast(int, LDBL_MIN_10_EXP); }
template <>        MC_TARGET_INLINE int mclimits_minexp10of<signed char>        (const signed char& x)        { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_minexp10of<short>              (const short& x)              { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_minexp10of<int>                (const int& x)                { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_minexp10of<long>               (const long& x)               { mc_unused(x); return 0;                             }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE int mclimits_minexp10of<long long>          (const long long& x)          { mc_unused(x); return 0;                             }
#	endif
template <>        MC_TARGET_INLINE int mclimits_minexp10of<unsigned char>      (const unsigned char& x)      { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_minexp10of<unsigned short>     (const unsigned short& x)     { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_minexp10of<unsigned int>       (const unsigned int& x)       { mc_unused(x); return 0;                             }
template <>        MC_TARGET_INLINE int mclimits_minexp10of<unsigned long>      (const unsigned long& x)      { mc_unused(x); return 0;                             }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE int mclimits_minexp10of<unsigned long long> (const unsigned long long& x) { mc_unused(x); return 0;                             }
#	endif
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS int mclimits_minexp10of (const float x)              { mc_unused(x); return mc_cast(int, FLT_MIN_10_EXP);  }
MC_TARGET_ALIAS int mclimits_minexp10of (const double x)             { mc_unused(x); return mc_cast(int, DBL_MIN_10_EXP);  }
MC_TARGET_ALIAS int mclimits_minexp10of (const long double x)        { mc_unused(x); return mc_cast(int, LDBL_MIN_10_EXP); }
MC_TARGET_ALIAS int mclimits_minexp10of (const signed char x)        { mc_unused(x); return 0;                             }
MC_TARGET_ALIAS int mclimits_minexp10of (const short x)              { mc_unused(x); return 0;                             }
MC_TARGET_ALIAS int mclimits_minexp10of (const int x)                { mc_unused(x); return 0;                             }
MC_TARGET_ALIAS int mclimits_minexp10of (const long x)               { mc_unused(x); return 0;                             }
#	if MC_TARGET_C99
MC_TARGET_ALIAS int mclimits_minexp10of (const long long x)          { mc_unused(x); return 0;                             }
#	endif
MC_TARGET_ALIAS int mclimits_minexp10of (const unsigned char x)      { mc_unused(x); return 0;                             }
MC_TARGET_ALIAS int mclimits_minexp10of (const unsigned short x)     { mc_unused(x); return 0;                             }
MC_TARGET_ALIAS int mclimits_minexp10of (const unsigned int x)       { mc_unused(x); return 0;                             }
MC_TARGET_ALIAS int mclimits_minexp10of (const unsigned long x)      { mc_unused(x); return 0;                             }
#	if MC_TARGET_C99
MC_TARGET_ALIAS int mclimits_minexp10of (const unsigned long long x) { mc_unused(x); return 0;                             }
#	endif
#	elif MC_TARGET_C11
#	define mclimits_minexp10of(x) mc_cast(int, _Generic(x \
	, float              : FLT_MIN_10_EXP                 \
	, double             : DBL_MIN_10_EXP                 \
	, long double        : LDBL_MIN_10_EXP                \
	, signed char        : 0                              \
	, short              : 0                              \
	, int                : 0                              \
	, long               : 0                              \
	, long long          : 0                              \
	, unsigned char      : 0                              \
	, unsigned short     : 0                              \
	, unsigned int       : 0                              \
	, unsigned long      : 0                              \
	, unsigned long long : 0                              \
))
#	elif MC_TARGET_HAVE_TYPEOF
#	if MC_TARGET_C99
#	define mclimits_minexp10of(x) mc_cast(int,                                         \
	(                                                                                  \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? FLT_MIN_10_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? DBL_MIN_10_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? LDBL_MIN_10_EXP \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long long)          ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long long) ? 0               \
		: 0                                                                             \
	))
#	else
#	define mclimits_minexp10of(x) mc_cast(int,                                         \
	(                                                                                  \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? FLT_MIN_10_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? DBL_MIN_10_EXP  \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? LDBL_MIN_10_EXP \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? 0               \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? 0               \
		: 0                                                                             \
	))
#	endif
#	else
#	define mclimits_minexp10of(x) (0)
#	endif

#pragma mark - mclimits_digits -

#	if MC_TARGET_CPP98
template <class T> MC_TARGET_INLINE int mclimits_digits                     (const T& x)                  { mc_unused(x); return 0;                                                        }
template <>        MC_TARGET_INLINE int mclimits_digits<float>              (const float& x)              { mc_unused(x); return mc_cast_expr(int, FLT_MANT_DIG);                          }
template <>        MC_TARGET_INLINE int mclimits_digits<double>             (const double& x)             { mc_unused(x); return mc_cast_expr(int, DBL_MANT_DIG);                          }
template <>        MC_TARGET_INLINE int mclimits_digits<long double>        (const long double& x)        { mc_unused(x); return mc_cast_expr(int, LDBL_MANT_DIG);                         }
template <>        MC_TARGET_INLINE int mclimits_digits<signed char>        (const signed char& x)        { mc_unused(x); return mc_cast_expr(int, CHAR_BIT - 1);                          }
template <>        MC_TARGET_INLINE int mclimits_digits<short>              (const short& x)              { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(short) - 1);          }
template <>        MC_TARGET_INLINE int mclimits_digits<int>                (const int& x)                { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(int)   - 1);          }
template <>        MC_TARGET_INLINE int mclimits_digits<long>               (const long& x)               { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(long)  - 1);          }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE int mclimits_digits<long long>          (const long long& x)          { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(long long) - 1);      }
#	endif
template <>        MC_TARGET_INLINE int mclimits_digits<unsigned char>      (const unsigned char& x)      { mc_unused(x); return mc_cast_expr(int, CHAR_BIT);                              }
template <>        MC_TARGET_INLINE int mclimits_digits<unsigned short>     (const unsigned short& x)     { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(unsigned short));     }
template <>        MC_TARGET_INLINE int mclimits_digits<unsigned int>       (const unsigned int& x)       { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(unsigned int));       }
template <>        MC_TARGET_INLINE int mclimits_digits<unsigned long>      (const unsigned long& x)      { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(unsigned long));      }
#	if MC_TARGET_CPP11
template <>        MC_TARGET_INLINE int mclimits_digits<unsigned long long> (const unsigned long long& x) { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(unsigned long long)); }
#	endif
#	elif MC_TARGET_HAVE_OVERLOADABLE
MC_TARGET_ALIAS int mclimits_digits (const float x)              { mc_unused(x); return mc_cast_expr(int, FLT_MANT_DIG);                          }
MC_TARGET_ALIAS int mclimits_digits (const double x)             { mc_unused(x); return mc_cast_expr(int, DBL_MANT_DIG);                          }
MC_TARGET_ALIAS int mclimits_digits (const long double x)        { mc_unused(x); return mc_cast_expr(int, LDBL_MANT_DIG);                         }
MC_TARGET_ALIAS int mclimits_digits (const signed char x)        { mc_unused(x); return mc_cast_expr(int, CHAR_BIT - 1);                          }
MC_TARGET_ALIAS int mclimits_digits (const short x)              { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(short) - 1);          }
MC_TARGET_ALIAS int mclimits_digits (const int x)                { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(int)   - 1);          }
MC_TARGET_ALIAS int mclimits_digits (const long x)               { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(long)  - 1);          }
#	if MC_TARGET_C99
MC_TARGET_ALIAS int mclimits_digits (const long long x)          { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(long long) - 1);      }
#	endif
MC_TARGET_ALIAS int mclimits_digits (const unsigned char x)      { mc_unused(x); return mc_cast_expr(int, CHAR_BIT);                              }
MC_TARGET_ALIAS int mclimits_digits (const unsigned short x)     { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(unsigned short));     }
MC_TARGET_ALIAS int mclimits_digits (const unsigned int x)       { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(unsigned int));       }
MC_TARGET_ALIAS int mclimits_digits (const unsigned long x)      { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(unsigned long));      }
#	if MC_TARGET_C99
MC_TARGET_ALIAS int mclimits_digits (const unsigned long long x) { mc_unused(x); return mc_cast_expr(int, CHAR_BIT * sizeof(unsigned long long)); }
#	endif
#	elif MC_TARGET_C11
#	define mclimits_digits(x) mc_cast(int, _Generic(x              \
	, float              : FLT_MANT_DIG                            \
	, double             : DBL_MANT_DIG                            \
	, long double        : LDBL_MANT_DIG                           \
	, signed char        : (CHAR_BIT - 1)                          \
	, short              : (CHAR_BIT * sizeof(short)     - 1)      \
	, int                : (CHAR_BIT * sizeof(int)       - 1)      \
	, long               : (CHAR_BIT * sizeof(long)      - 1)      \
	, long long          : (CHAR_BIT * sizeof(long long) - 1)      \
	, unsigned char      : CHAR_BIT                                \
	, unsigned short     : (CHAR_BIT * sizeof(unsigned short))     \
	, unsigned int       : (CHAR_BIT * sizeof(unsigned int))       \
	, unsigned long      : (CHAR_BIT * sizeof(unsigned long))      \
	, unsigned long long : (CHAR_BIT * sizeof(unsigned long long)) \
))
#	elif MC_TARGET_HAVE_TYPEOF
#	if MC_TARGET_C99
#	define mclimits_digits(x) mc_cast(int,                                                                     \
	(                                                                                                          \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? FLT_MANT_DIG                            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? DBL_MANT_DIG                            \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? LDBL_MANT_DIG)                          \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? (CHAR_BIT - 1)                          \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? (CHAR_BIT * sizeof(short)     - 1)      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? (CHAR_BIT * sizeof(int)       - 1)      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? (CHAR_BIT * sizeof(long)      - 1)      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long long)          ? (CHAR_BIT * sizeof(long long) - 1)      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? CHAR_BIT                                \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? (CHAR_BIT * sizeof(unsigned short))     \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? (CHAR_BIT * sizeof(unsigned int))       \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? (CHAR_BIT * sizeof(unsigned long))      \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long long) ? (CHAR_BIT * sizeof(unsigned long long)) \
		: 0                                                                                                     \
	))
#	else
#	define mclimits_digits(x) mc_cast(int,                                                                   \
	(                                                                                                        \
		  MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), float)              ? FLT_MANT_DIG                          \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), double)             ? DBL_MANT_DIG                          \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long double)        ? LDBL_MANT_DIG                         \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), signed char)        ? (CHAR_BIT - 1)                        \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), short)              ? (CHAR_BIT * sizeof(short)     - 1)    \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), int)                ? (CHAR_BIT * sizeof(int)       - 1)    \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), long)               ? (CHAR_BIT * sizeof(long)      - 1)    \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned char)      ? CHAR_BIT                              \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned short)     ? (CHAR_BIT * sizeof(unsigned short))   \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned int)       ? (CHAR_BIT * sizeof(unsigned int))     \
		: MC_TARGET_TYPEISOF(MC_TARGET_TYPEOF(x), unsigned long)      ? (CHAR_BIT * sizeof(unsigned long))    \
		: 0                                                                                                   \
	))
#	endif
#	else
#	define mclimits_digits(x) (0)
#	endif

#endif /* !MCLIMITS_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_sincospi.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_remint2.h>
#include <macadam/details/math/mc_sincos.h>

#ifndef MC_SINCOSPI_H
#define MC_SINCOSPI_H

#pragma mark - mc_sincospi -

MC_TARGET_FUNC void mc_sincospif(const float x, float * sinp, float * cosp)
{
	float z   = 0.0f, ss = 0.0f, cc = 0.0f;
	int64_t i = mc_remint2f(x, &z) & 3;
	z         = MCK_KF(MCK_PI) * z;
	mc_sincosf(z, &ss, &cc);
	switch (i)
	{
		case 0:
			*sinp = ss;
			*cosp = cc;
		break;
		case 1:
			*sinp =  cc;
			*cosp = -ss;
		break;
		case 2:
			*sinp = -ss;
			*cosp = -cc;
		break;
		default:
			*sinp = -cc;
			*cosp =  ss;
	}
}

MC_TARGET_FUNC void mc_sincospi(const double x, double * sinp, double * cosp)
{
	double z  = 0.0, ss = 0.0, cc = 0.0;
	int64_t i = mc_remint2(x, &z) & 3;
	z         = MCK_K(MCK_PI) * z;
	mc_sincos(z, &ss, &cc);
	switch (i)
	{
		case 0:
			*sinp = ss;
			*cosp = cc;
		break;
		case 1:
			*sinp =  cc;
			*cosp = -ss;
		break;
		case 2:
			*sinp = -ss;
			*cosp = -cc;
		break;
		default:
			*sinp = -cc;
			*cosp =  ss;
	}
}

MC_TARGET_FUNC void mc_sincospil(const long double x, long double * sinp, long double * cosp)
{
	long double z = 0.0L, ss = 0.0L, cc = 0.0L;
	int64_t i     = mc_remint2l(x, &z) & 3;
#	if (MC_TARGET_C99 || MC_TARGET_CPP17) && defined(M_PIl)
	z             = M_PIl * z;
#	else
	z             = MCK_KL(MCK_PI) * z;
#	endif
	mc_sincosl(z, &ss, &cc);
	switch (i)
	{
		case 0:
			*sinp = ss;
			*cosp = cc;
		break;
		case 1:
			*sinp =  cc;
			*cosp = -ss;
		break;
		case 2:
			*sinp = -ss;
			*cosp = -cc;
		break;
		default:
			*sinp = -cc;
			*cosp =  ss;
	}
}

#endif /* !MC_SINCOSPI_H */

/* EOF */
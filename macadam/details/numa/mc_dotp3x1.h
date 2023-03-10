//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_dotp3x1.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_DOTP3X1_H
#define MC_DOTP3X1_H

#pragma mark - mc_dotp3x1 -

MC_TARGET_FUNC float mc_dotp3x1f(const int n, const int j, const int k, const float * a, const float * b, const int f)
{
	float s = 0.0f, e = 0.0f, x, y;

	switch (f)
	{
		case 0:
			s = s + (a[j] * b[k]);
			s = s + (a[n + j] * b[n + k]);
			s = s + (a[(n * 2) + j] * b[(n * 2) + k]);
		break;
		case 1:
			mc_twoproductf(a[j], b[k], &x, &y);
			mc_twosumf(s, x + y, &s, &y);
			e = e + y;

			mc_twoproductf(a[n + j], b[n + k], &x, &y);
			mc_twosumf(s, x + y, &s, &y);
			e = e + y;

			mc_twoproductf(a[(n * 2) + j], b[(n * 2) + k], &x, &y);
			mc_twosumf(s, x + y, &s, &y);
			e = e + y;
		break;
	}
	return s + e;
}

MC_TARGET_FUNC double mc_dotp3x1ff(const int n, const int j, const int k, const float * a, const float * b, const int f)
{
	double s = 0.0, e = 0.0, x, y;

	switch (f)
	{
		case 0:
			s = s + (mc_cast(double, a[j]) * mc_cast(double, b[k]));
			s = s + (mc_cast(double, a[n + j]) * mc_cast(double, b[n + k]));
			s = s + (mc_cast(double, a[(n * 2) + j]) * mc_cast(double, b[(n * 2) + k]));
		break;
		case 1:
			mc_twoproduct(mc_cast(double, a[j]), mc_cast(double, b[k]), &x, &y);
			mc_twosum(s, x + y, &s, &y);
			e = e + y;

			mc_twoproduct(mc_cast(double, a[n + j]), mc_cast(double, b[n + k]), &x, &y);
			mc_twosum(s, x + y, &s, &y);
			e = e + y;

			mc_twoproduct(mc_cast(double, a[(n * 2) + j]), mc_cast(double, b[(n * 2) + k]), &x, &y);
			mc_twosum(s, x + y, &s, &y);
			e = e + y;
		break;
	}
	return s + e;
}

MC_TARGET_FUNC double mc_dotp3x1(const int n, const int j, const int k, const double * a, const double * b, const int f)
{
	double s = 0.0, e = 0.0, x, y;

	switch (f)
	{
		case 0:
			s = s + (a[j] * b[k]);
			s = s + (a[n + j] * b[n + k]);
			s = s + (a[(n * 2) + j] * b[(n * 2) + k]);
		break;
		case 1:
			mc_twoproduct(a[j], b[k], &x, &y);
			mc_twosum(s, x + y, &s, &y);
			e = e + y;

			mc_twoproduct(a[n + j], b[n + k], &x, &y);
			mc_twosum(s, x + y, &s, &y);
			e = e + y;

			mc_twoproduct(a[(n * 2) + j], b[(n * 2) + k], &x, &y);
			mc_twosum(s, x + y, &s, &y);
			e = e + y;
		break;
	}
	return s + e;
}

MC_TARGET_FUNC long double mc_dotp3x1l(const int n, const int j, const int k, const long double * a, const long double * b, const int f)
{
	long double s = 0.0L, e = 0.0L, x, y;

	switch (f)
	{
		case 0:
			s = s + (a[j] * b[k]);
			s = s + (a[n + j] * b[n + k]);
			s = s + (a[(n * 2) + j] * b[(n * 2) + k]);
		break;
		case 1:
			mc_twoproductl(a[j], b[k], &x, &y);
			mc_twosuml(s, x + y, &s, &y);
			e = e + y;

			mc_twoproductl(a[n + j], b[n + k], &x, &y);
			mc_twosuml(s, x + y, &s, &y);
			e = e + y;

			mc_twoproductl(a[(n * 2) + j], b[(n * 2) + k], &x, &y);
			mc_twosuml(s, x + y, &s, &y);
			e = e + y;
		break;
	}
	return s + e;
}

#endif /* !MC_DOTP3X1_H */

/* EOF */
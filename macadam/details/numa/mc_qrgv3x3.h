//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_qrgv3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/math/mc_hypot2.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_fmax.h>

#ifndef MC_QRGV3X3_H
#define MC_QRGV3X3_H

#pragma mark - mc_gvrot -

MC_TARGET_FUNC void mc_gvrotf(float a1, float a2, float tol, float * ch, float * sh, float * r)
{
//!# Givens rotation.
	const int wantr = mc_nonnullptr(r);
	float w;

//!# Sanity check.
	 w  = mc_hypot2f(a1, a2);
	*sh = w > tol ? a2 : 0.0f;
	*ch = mc_fabsf(a1) + mc_fmaxf(w, tol);

//!# Sign check.
	 w  = *sh;
	*sh = mc_copysignf(1.0f, a1) < 0.0f ? *ch : *sh;
	*ch = mc_copysignf(1.0f, a1) < 0.0f ?  w  : *ch;

//!# Assigning r without normalization.
	if (wantr) {
		*r = mc_hypot2f(*ch, *sh);
	} else {
//!# Normalizing ch and sh by r.
		w = mc_hypot2f(*ch, *sh);
		if (w != 0.0f) {
			 w  =  1.0f / w;
			*ch = *ch * w;
			*sh = *sh * w;
		}
	}
}

MC_TARGET_FUNC void mc_gvrotff(float a1, float a2, float tol, double * ch, double * sh, double * r)
{
//!# Givens rotation.
	const int wantr = mc_nonnullptr(r);
	double w, a1d, a2d, told;

	a1d  = mc_cast(double, a1);
	a2d  = mc_cast(double, a2);
	told = mc_cast(double, tol);

//!# Sanity check.
	 w   = mc_hypot2(a1d, a2d);
	*sh  = w > told ? a2d : 0.0;
	*ch  = mc_fabs(a1d) + mc_fmax(w, told);

//!# Sign check.
	 w   = *sh;
	*sh  = mc_copysign(1.0, a1d) < 0.0 ? *ch : *sh;
	*ch  = mc_copysign(1.0, a1d) < 0.0 ?  w  : *ch;

//!# Assigning r without normalization.
	if (wantr) {
		*r = mc_hypot2(*ch, *sh);
	} else {
//!# Normalizing ch and sh by r.
		w = mc_hypot2(*ch, *sh);
		if (w != 0.0) {
			 w  =  1.0 / w;
			*ch = *ch * w;
			*sh = *sh * w;
		}
	}
}

MC_TARGET_FUNC void mc_gvrot(double a1, double a2, double tol, double * ch, double * sh, double * r)
{
//!# Givens rotation.
	const int wantr = mc_nonnullptr(r);
	double w;

//!# Sanity check.
	 w  = mc_hypot2(a1, a2);
	*sh = w > tol ? a2 : 0.0;
	*ch = mc_fabs(a1) + mc_fmax(w, tol);

//!# Sign check.
	 w  = *sh;
	*sh = mc_copysign(1.0, a1) < 0.0 ? *ch : *sh;
	*ch = mc_copysign(1.0, a1) < 0.0 ?  w  : *ch;

//!# Assigning r without normalization.
	if (wantr) {
		*r = mc_hypot2(*ch, *sh);
	} else {
//!# Normalizing ch and sh by r.
		w = mc_hypot2(*ch, *sh);
		if (w != 0.0) {
			 w  =  1.0 / w;
			*ch = *ch * w;
			*sh = *sh * w;
		}
	}
}

MC_TARGET_FUNC void mc_gvrotl(long double a1, long double a2, long double tol, long double * ch, long double * sh, long double * r)
{
//!# Givens rotation.
	const int wantr = mc_nonnullptr(r);
	long double w;

//!# Sanity check.
	 w  = mc_hypot2l(a1, a2);
	*sh = w > tol ? a2 : 0.0L;
	*ch = mc_fabsl(a1) + mc_fmaxl(w, tol);

//!# Sign check.
	 w  = *sh;
	*sh = mc_copysignl(1.0L, a1) < 0.0L ? *ch : *sh;
	*ch = mc_copysignl(1.0L, a1) < 0.0L ?  w  : *ch;

//!# Assigning r without normalization.
	if (wantr) {
		*r = mc_hypot2l(*ch, *sh);
	} else {
//!# Normalizing ch and sh by r.
		w = mc_hypot2l(*ch, *sh);
		if (w != 0.0L) {
			 w  =  1.0L / w;
			*ch = *ch * w;
			*sh = *sh * w;
		}
	}
}

#pragma mark - mc_qrgv3x3 -

MC_TARGET_FUNC int mc_qrgv3x3f(const float a[9], float q[9], float r[9])
{
//!# A and Q may be the same. Using Givens rotations method.
	const float tol = 2.0f * MCLIMITS_EPSILONF;

	float a0, b0, ch1, sh1, ch2, sh2, ch3, sh3, sh12, sh22, sh32;

	float a11 = a[0], a12 = a[1], a13 = a[2];
	float a21 = a[3], a22 = a[4], a23 = a[5];
	float a31 = a[6], a32 = a[7], a33 = a[8];

//!# First givens rotation (ch, 0, 0, sh)
	mc_gvrotf(a11, a21, tol, &ch1, &sh1, MC_NULLPTR);
	a0 = 1.0f - 2.0f * sh1 * sh1;
	b0 = 2.0f        * ch1 * sh1;

//!# Computing B=Q'*B
	r[0] =  a0 * a11 + b0 * a21; r[1] =  a0 * a12 + b0 * a22; r[2] =  a0 * a13 + b0 * a23;
	r[3] = -b0 * a11 + a0 * a21; r[4] = -b0 * a12 + a0 * a22; r[5] = -b0 * a13 + a0 * a23;
	r[6] =  a31;                 r[7] =  a32;                 r[8] =  a33;

//!# Second givens rotation (ch, 0, -sh, 0)
	mc_gvrotf(r[0], r[6], tol, &ch2, &sh2, MC_NULLPTR);
	a0 = 1.0f - 2.0f * sh2 * sh2;
	b0 = 2.0f        * ch2 * sh2;

//!# Computing B=Q'*B
	a11 =  a0 * r[0] + b0 * r[6]; a12 =  a0 * r[1] + b0 * r[7]; a13 =  a0 * r[2] + b0 * r[8];
	a21 =  r[3];                  a22 =  r[4];                  a23 =  r[5];
	a31 = -b0 * r[0] + a0 * r[6]; a32 = -b0 * r[1] + a0 * r[7]; a33 = -b0 * r[2] + a0 * r[8];

//!# Third givens rotation (ch, sh, 0, 0)
	mc_gvrotf(a22, a32, tol, &ch3, &sh3, MC_NULLPTR);
	a0 = 1.0f - 2.0f * sh3 * sh3;
	b0 = 2.0f        * ch3 * sh3;

//!# Computing and finalizing R.
	r[0] =  a11;                 r[1] =  a12;                 r[2] =  a13;
	r[3] =  a0 * a21 + b0 * a31; r[4] =  a0 * a22 + b0 * a32; r[5] =  a0 * a23 + b0 * a33;
	r[6] = -b0 * a21 + a0 * a31; r[7] = -b0 * a22 + a0 * a32; r[8] = -b0 * a23 + a0 * a33;

	sh12 = sh1 * sh1;
	sh22 = sh2 * sh2;
	sh32 = sh3 * sh3;

//!# Computing Q such as Q=Q1*Q2*Q3
	q[0] =  (-1.0f + 2.0f * sh12) * (-1.0f + 2.0f * sh22); 
	q[1] =  4.0f * ch2 * ch3 * (-1.0f + 2.0f * sh12) * sh2 * sh3 + 2.0f * ch1 * sh1 * (-1.0f + 2.0f * sh32); 
	q[2] =  4.0f * ch1 * ch3 * sh1 * sh3 - 2.0f * ch2 * (-1.0f + 2.0f * sh12) * sh2 * (-1.0f + 2.0f * sh32);

	q[3] =  2.0f * ch1 * sh1 * (1.0f - 2.0f * sh22); 
	q[4] = -8.0f * ch1 * ch2 * ch3 * sh1 * sh2 * sh3 + (-1.0f + 2.0f * sh12) * (-1.0f + 2.0f * sh32); 
	q[5] = -2.0f * ch3 * sh3 + 4.0f * sh1 * (ch3 * sh1 * sh3+ch1 * ch2 * sh2 * (-1.0f + 2.0f * sh32));

	q[6] =  2.0f * ch2 * sh2; 
	q[7] =  2.0f * ch3 * (1.0f - 2.0f * sh22) * sh3; 
	q[8] =  (-1.0f + 2.0f * sh22) * (-1.0f + 2.0f * sh32);

	return 0;
}

MC_TARGET_PROC int mc_qrgv3x3ff(const float a[9], double q[9], double r[9])
{
//!# Using Givens rotations method.
	const double tol = 2.0 * MCLIMITS_EPSILON;

	double a0, b0, ch1, sh1, ch2, sh2, ch3, sh3, sh12, sh22, sh32;

	double a11 = mc_cast(double, a[0]), a12 = mc_cast(double, a[1]), a13 = mc_cast(double, a[2]);
	double a21 = mc_cast(double, a[3]), a22 = mc_cast(double, a[4]), a23 = mc_cast(double, a[5]);
	double a31 = mc_cast(double, a[6]), a32 = mc_cast(double, a[7]), a33 = mc_cast(double, a[8]);

//!# First givens rotation (ch, 0, 0, sh)
	mc_gvrot(a11, a21, tol, &ch1, &sh1, MC_NULLPTR);
	a0 = 1.0 - 2.0 * sh1 * sh1;
	b0 = 2.0       * ch1 * sh1;

//!# Computing B=Q'*B
	r[0] =  a0 * a11 + b0 * a21; r[1] =  a0 * a12 + b0 * a22; r[2] =  a0 * a13 + b0 * a23;
	r[3] = -b0 * a11 + a0 * a21; r[4] = -b0 * a12 + a0 * a22; r[5] = -b0 * a13 + a0 * a23;
	r[6] =  a31;                 r[7] =  a32;                 r[8] =  a33;

//!# Second givens rotation (ch, 0, -sh, 0)
	mc_gvrot(r[0], r[6], tol, &ch2, &sh2, MC_NULLPTR);
	a0 = 1.0 - 2.0 * sh2 * sh2;
	b0 = 2.0       * ch2 * sh2;

//!# Computing B=Q'*B
	a11 =  a0 * r[0] + b0 * r[6]; a12 =  a0 * r[1] + b0 * r[7]; a13 =  a0 * r[2] + b0 * r[8];
	a21 =  r[3];                  a22 =  r[4];                  a23 =  r[5];
	a31 = -b0 * r[0] + a0 * r[6]; a32 = -b0 * r[1] + a0 * r[7]; a33 = -b0 * r[2] + a0 * r[8];

//!# Third givens rotation (ch, sh, 0, 0)
	mc_gvrot(a22, a32, tol, &ch3, &sh3, MC_NULLPTR);
	a0 = 1.0 - 2.0 * sh3 * sh3;
	b0 = 2.0       * ch3 * sh3;

//!# Computing and finalizing R.
	r[0] =  a11;                 r[1] =  a12;                 r[2] =  a13;
	r[3] =  a0 * a21 + b0 * a31; r[4] =  a0 * a22 + b0 * a32; r[5] =  a0 * a23 + b0 * a33;
	r[6] = -b0 * a21 + a0 * a31; r[7] = -b0 * a22 + a0 * a32; r[8] = -b0 * a23 + a0 * a33;

	sh12 = sh1 * sh1;
	sh22 = sh2 * sh2;
	sh32 = sh3 * sh3;

//!# Computing Q such as Q=Q1*Q2*Q3
	q[0] =  (-1.0 + 2.0 * sh12) * (-1.0 + 2.0 * sh22); 
	q[1] =  4.0 * ch2 * ch3 * (-1.0 + 2.0 * sh12) * sh2 * sh3 + 2.0 * ch1 * sh1 * (-1.0 + 2.0 * sh32); 
	q[2] =  4.0 * ch1 * ch3 * sh1 * sh3 - 2.0 * ch2 * (-1.0 + 2.0 * sh12) * sh2 * (-1.0 + 2.0 * sh32);

	q[3] =  2.0 * ch1 * sh1 * (1.0 - 2.0 * sh22); 
	q[4] = -8.0 * ch1 * ch2 * ch3 * sh1 * sh2 * sh3 + (-1.0 + 2.0 * sh12) * (-1.0 + 2.0 * sh32); 
	q[5] = -2.0 * ch3 * sh3 + 4.0 * sh1 * (ch3 * sh1 * sh3+ch1 * ch2 * sh2 * (-1.0 + 2.0 * sh32));

	q[6] =  2.0 * ch2 * sh2; 
	q[7] =  2.0 * ch3 * (1.0 - 2.0 * sh22) * sh3; 
	q[8] =  (-1.0 + 2.0 * sh22) * (-1.0 + 2.0 * sh32);

	return 0;
}

MC_TARGET_PROC int mc_qrgv3x3(const double a[9], double q[9], double r[9])
{
//!# A and Q may be the same. Using Givens rotations method.
	const double tol = 2.0 * MCLIMITS_EPSILON;

	double a0, b0, ch1, sh1, ch2, sh2, ch3, sh3, sh12, sh22, sh32;

	double a11 = a[0], a12 = a[1], a13 = a[2];
	double a21 = a[3], a22 = a[4], a23 = a[5];
	double a31 = a[6], a32 = a[7], a33 = a[8];

//!# First givens rotation (ch, 0, 0, sh)
	mc_gvrot(a11, a21, tol, &ch1, &sh1, MC_NULLPTR);
	a0 = 1.0 - 2.0 * sh1 * sh1;
	b0 = 2.0       * ch1 * sh1;

//!# Computing B=Q'*B
	r[0] =  a0 * a11 + b0 * a21; r[1] =  a0 * a12 + b0 * a22; r[2] =  a0 * a13 + b0 * a23;
	r[3] = -b0 * a11 + a0 * a21; r[4] = -b0 * a12 + a0 * a22; r[5] = -b0 * a13 + a0 * a23;
	r[6] =  a31;                 r[7] =  a32;                 r[8] =  a33;

//!# Second givens rotation (ch, 0, -sh, 0)
	mc_gvrot(r[0], r[6], tol, &ch2, &sh2, MC_NULLPTR);
	a0 = 1.0 - 2.0 * sh2 * sh2;
	b0 = 2.0       * ch2 * sh2;

//!# Computing B=Q'*B
	a11 =  a0 * r[0] + b0 * r[6]; a12 =  a0 * r[1] + b0 * r[7]; a13 =  a0 * r[2] + b0 * r[8];
	a21 =  r[3];                  a22 =  r[4];                  a23 =  r[5];
	a31 = -b0 * r[0] + a0 * r[6]; a32 = -b0 * r[1] + a0 * r[7]; a33 = -b0 * r[2] + a0 * r[8];

//!# Third givens rotation (ch, sh, 0, 0)
	mc_gvrot(a22, a32, tol, &ch3, &sh3, MC_NULLPTR);
	a0 = 1.0 - 2.0 * sh3 * sh3;
	b0 = 2.0       * ch3 * sh3;

//!# Computing and finalizing R.
	r[0] =  a11;                 r[1] =  a12;                 r[2] =  a13;
	r[3] =  a0 * a21 + b0 * a31; r[4] =  a0 * a22 + b0 * a32; r[5] =  a0 * a23 + b0 * a33;
	r[6] = -b0 * a21 + a0 * a31; r[7] = -b0 * a22 + a0 * a32; r[8] = -b0 * a23 + a0 * a33;

	sh12 = sh1 * sh1;
	sh22 = sh2 * sh2;
	sh32 = sh3 * sh3;

//!# Computing Q such as Q=Q1*Q2*Q3
	q[0] =  (-1.0 + 2.0 * sh12) * (-1.0 + 2.0 * sh22); 
	q[1] =  4.0 * ch2 * ch3 * (-1.0 + 2.0 * sh12) * sh2 * sh3 + 2.0 * ch1 * sh1 * (-1.0 + 2.0 * sh32); 
	q[2] =  4.0 * ch1 * ch3 * sh1 * sh3 - 2.0 * ch2 * (-1.0 + 2.0 * sh12) * sh2 * (-1.0 + 2.0 * sh32);

	q[3] =  2.0 * ch1 * sh1 * (1.0 - 2.0 * sh22); 
	q[4] = -8.0 * ch1 * ch2 * ch3 * sh1 * sh2 * sh3 + (-1.0 + 2.0 * sh12) * (-1.0 + 2.0 * sh32); 
	q[5] = -2.0 * ch3 * sh3 + 4.0 * sh1 * (ch3 * sh1 * sh3+ch1 * ch2 * sh2 * (-1.0 + 2.0 * sh32));

	q[6] =  2.0 * ch2 * sh2; 
	q[7] =  2.0 * ch3 * (1.0 - 2.0 * sh22) * sh3; 
	q[8] =  (-1.0 + 2.0 * sh22) * (-1.0 + 2.0 * sh32);

	return 0;
}

MC_TARGET_PROC int mc_qrgv3x3l(const long double a[9], long double q[9], long double r[9])
{
//!# A and Q may be the same. Using Givens rotations method.
	const long double tol = 2.0L * MCLIMITS_EPSILONL;

	long double a0, b0, ch1, sh1, ch2, sh2, ch3, sh3, sh12, sh22, sh32;

	long double a11 = a[0], a12 = a[1], a13 = a[2];
	long double a21 = a[3], a22 = a[4], a23 = a[5];
	long double a31 = a[6], a32 = a[7], a33 = a[8];

//!# First givens rotation (ch, 0, 0, sh)
	mc_gvrotl(a11, a21, tol, &ch1, &sh1, MC_NULLPTR);
	a0 = 1.0L - 2.0L * sh1 * sh1;
	b0 = 2.0L        * ch1 * sh1;

//!# Computing B=Q'*B
	r[0] =  a0 * a11 + b0 * a21; r[1] =  a0 * a12 + b0 * a22; r[2] =  a0 * a13 + b0 * a23;
	r[3] = -b0 * a11 + a0 * a21; r[4] = -b0 * a12 + a0 * a22; r[5] = -b0 * a13 + a0 * a23;
	r[6] =  a31;                 r[7] =  a32;                 r[8] =  a33;

//!# Second givens rotation (ch, 0, -sh, 0)
	mc_gvrotl(r[0], r[6], tol, &ch2, &sh2, MC_NULLPTR);
	a0 = 1.0L - 2.0L * sh2 * sh2;
	b0 = 2.0L        * ch2 * sh2;

//!# Computing B=Q'*B
	a11 =  a0 * r[0] + b0 * r[6]; a12 =  a0 * r[1] + b0 * r[7]; a13 =  a0 * r[2] + b0 * r[8];
	a21 =  r[3];                  a22 =  r[4];                  a23 =  r[5];
	a31 = -b0 * r[0] + a0 * r[6]; a32 = -b0 * r[1] + a0 * r[7]; a33 = -b0 * r[2] + a0 * r[8];

//!# Third givens rotation (ch, sh, 0, 0)
	mc_gvrotl(a22, a32, tol, &ch3, &sh3, MC_NULLPTR);
	a0 = 1.0L - 2.0L * sh3 * sh3;
	b0 = 2.0L        * ch3 * sh3;

//!# Computing and finalizing R.
	r[0] =  a11;                 r[1] =  a12;                 r[2] =  a13;
	r[3] =  a0 * a21 + b0 * a31; r[4] =  a0 * a22 + b0 * a32; r[5] =  a0 * a23 + b0 * a33;
	r[6] = -b0 * a21 + a0 * a31; r[7] = -b0 * a22 + a0 * a32; r[8] = -b0 * a23 + a0 * a33;

	sh12 = sh1 * sh1;
	sh22 = sh2 * sh2;
	sh32 = sh3 * sh3;

//!# Computing Q such as Q=Q1*Q2*Q3
	q[0] =  (-1.0L + 2.0L * sh12) * (-1.0L + 2.0L * sh22); 
	q[1] =  4.0L * ch2 * ch3 * (-1.0L + 2.0L * sh12) * sh2 * sh3 + 2.0L * ch1 * sh1 * (-1.0L + 2.0L * sh32); 
	q[2] =  4.0L * ch1 * ch3 * sh1 * sh3 - 2.0L * ch2 * (-1.0L + 2.0L * sh12) * sh2 * (-1.0L + 2.0L * sh32);

	q[3] =  2.0L * ch1 * sh1 * (1.0L - 2.0L * sh22); 
	q[4] = -8.0L * ch1 * ch2 * ch3 * sh1 * sh2 * sh3 + (-1.0L + 2.0L * sh12) * (-1.0L + 2.0L * sh32); 
	q[5] = -2.0L * ch3 * sh3 + 4.0L * sh1 * (ch3 * sh1 * sh3+ch1 * ch2 * sh2 * (-1.0L + 2.0L * sh32));

	q[6] =  2.0L * ch2 * sh2; 
	q[7] =  2.0L * ch3 * (1.0L - 2.0L * sh22) * sh3; 
	q[8] =  (-1.0L + 2.0L * sh22) * (-1.0L + 2.0L * sh32);

	return 0;
}

#endif /* !MC_QRGV3X3_H */

/* EOF */
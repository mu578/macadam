//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lup3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/numa/mc_eye3x3.h>
#include <macadam/mcswap.h>

#ifndef MC_LUP3X3_H
#define MC_LUP3X3_H

#pragma mark - mc_lup3x3 -

MC_TARGET_FUNC int mc_lup3x3f(const float a[9], float l[9], float u[9], float p[9], int pvi[3])
{
//!# A and L may be the same. Using a closed-form expression.
//!# Returns A=L(DU) as per Doolittle's method.
	const int wantpvi = mc_nonnullptr(pvi);
	const int wantp   = mc_nonnullptr(p);

	float w;
	int p0 = 0, p1 = 1, p2 = 2, j, pv = 0;

	float a11, a12, a13;
	float a21, a22, a23;
	float a31, a32, a33;

	float e111 =  1.0f, e112 = 0.0f, e113 = 0.0f;
	float e121 = -0.0f, e122 = 1.0f, e123 = 0.0f;
	float e131 = -0.0f, e132 = 0.0f, e133 = 1.0f;

	float e211 = 1.0f, e212 =  0.0f, e213 = 0.0f;
	float e221 = 0.0f, e222 =  1.0f, e223 = 0.0f;
	float e231 = 0.0f, e232 = -0.0f, e233 = 1.0f;

	float e1a11, e1a12, e1a13;
	float e1a21, e1a22, e1a23;
	float e1a31, e1a32, e1a33;

//!# Copying a.
	a11 = a[0]; a12 = a[1]; a13 = a[2];
	a21 = a[3]; a22 = a[4]; a23 = a[5];
	a31 = a[6]; a32 = a[7]; a33 = a[8];
//!# Finding largest pivot bottom to top.
	if (mc_fabsf(a11) < mc_fabsf(a31) && mc_fabsf(a31) > mc_fabsf(a21)) {
		mcswap_var(j, p0, p2);
		mcswap_var(w, a11, a31);
		mcswap_var(w, a12, a32);
		mcswap_var(w, a13, a33);
		++pv;
	} else if (mc_fabsf(a11) < (a21)) {
		mcswap_var(j, p0, p1);
		mcswap_var(w, a11, a21);
		mcswap_var(w, a12, a22);
		mcswap_var(w, a13, a23);
		++pv;
	}
//!# If selected pivot value != 0.
	if (a11 != 0.0f) {
//!# Forming first elementary matrix.
		e121 = -(a21 / a11);
		e131 = -(a31 / a11);
	}
//!# Computing e1a=e1*a.
	e1a11 = (e111 * a11) + (e112 * a21) + (e113 * a31);
	e1a12 = (e111 * a12) + (e112 * a22) + (e113 * a32);
	e1a13 = (e111 * a13) + (e112 * a23) + (e113 * a33);

	e1a21 = (e121 * a11) + (e122 * a21) + (e123 * a31);
	e1a22 = (e121 * a12) + (e122 * a22) + (e123 * a32);
	e1a23 = (e121 * a13) + (e122 * a23) + (e123 * a33);

	e1a31 = (e131 * a11) + (e132 * a21) + (e133 * a31);
	e1a32 = (e131 * a12) + (e132 * a22) + (e133 * a32);
	e1a33 = (e131 * a13) + (e132 * a23) + (e133 * a33);

//!# Finding largest pivot bottom to top.
	if (mc_fabsf(e1a32) > mc_fabsf(e1a22)) {
		mcswap_var(j, p1, p2);
		mcswap_var(w, e1a21, e1a31);
		mcswap_var(w, e1a22, e1a32);
		mcswap_var(w, e1a23, e1a33);
		mcswap_var(w, e121, e131);
		++pv;
	}
//!# If selected pivot value != 0.
	if (e1a22 != 0.0f) {
//!# Forming second elementary matrix e2.
		e232 = -(e1a32 / e1a22);
	}
//!# Computing U such as U=e2*e1a;
	u[0] = (e211 * e1a11) + (e212 * e1a21) + (e213 * e1a31);
	u[1] = (e211 * e1a12) + (e212 * e1a22) + (e213 * e1a32);
	u[2] = (e211 * e1a13) + (e212 * e1a23) + (e213 * e1a33);

	u[3] = (e221 * e1a11) + (e222 * e1a21) + (e223 * e1a31);
	u[4] = (e221 * e1a12) + (e222 * e1a22) + (e223 * e1a32);
	u[5] = (e221 * e1a13) + (e222 * e1a23) + (e223 * e1a33);

	u[6] = (e231 * e1a11) + (e232 * e1a21) + (e233 * e1a31);
	u[7] = (e231 * e1a12) + (e232 * e1a22) + (e233 * e1a32);
	u[8] = (e231 * e1a13) + (e232 * e1a23) + (e233 * e1a33);

//!# Computing L from e1 and e2 such as e1^-1*e2^-1.
	l[0] =  1.0f; l[1] =  0.0f; l[2] = 0.0f;
	l[3] = -e121; l[4] =  1.0f; l[5] = 0.0f;
	l[6] = -e131; l[7] = -e232; l[8] = 1.0f;

//!# Computing P or assigning pivot indeces.
	if (wantpvi) {
		pvi[0] = p0; pvi[1] = p1; pvi[2] = p2; 
	} else if (wantp) {
		p[0] = p0 == 0 ? 1.0f : 0.0f; p[1] = p0 == 1 ? 1.0f : 0.0f; p[2] = p0 == 2 ? 1.0f : 0.0f;
		p[3] = p1 == 0 ? 1.0f : 0.0f; p[4] = p1 == 1 ? 1.0f : 0.0f; p[5] = p1 == 2 ? 1.0f : 0.0f;
		p[6] = p2 == 0 ? 1.0f : 0.0f; p[7] = p2 == 1 ? 1.0f : 0.0f; p[8] = p2 == 2 ? 1.0f : 0.0f;
	}
	return pv;
}

MC_TARGET_FUNC int mc_lup3x3ff(const float a[9], double l[9], double u[9], double p[9], int pvi[3])
{
//!# A and L may be the same. Using a closed-form expression.
//!# Returns A=L(DU) as per Doolittle's method.
	const int wantpvi = mc_nonnullptr(pvi);
	const int wantp   = mc_nonnullptr(p);

	double w;
	int p0 = 0, p1 = 1, p2 = 2, j, pv = 0;

	double a11, a12, a13;
	double a21, a22, a23;
	double a31, a32, a33;

	double e111 =  1.0, e112 = 0.0, e113 = 0.0;
	double e121 = -0.0, e122 = 1.0, e123 = 0.0;
	double e131 = -0.0, e132 = 0.0, e133 = 1.0;

	double e211 = 1.0, e212 =  0.0, e213 = 0.0;
	double e221 = 0.0, e222 =  1.0, e223 = 0.0;
	double e231 = 0.0, e232 = -0.0, e233 = 1.0;

	double e1a11, e1a12, e1a13;
	double e1a21, e1a22, e1a23;
	double e1a31, e1a32, e1a33;

//!# Copying a.
	a11 = mc_cast(double, a[0]); a12 = mc_cast(double, a[1]); a13 = mc_cast(double, a[2]);
	a21 = mc_cast(double, a[3]); a22 = mc_cast(double, a[4]); a23 = mc_cast(double, a[5]);
	a31 = mc_cast(double, a[6]); a32 = mc_cast(double, a[7]); a33 = mc_cast(double, a[8]);
//!# Finding largest pivot bottom to top.
	if (mc_fabs(a11) < mc_fabs(a31) && mc_fabs(a31) > mc_fabs(a21)) {
		mcswap_var(j, p0, p2);
		mcswap_var(w, a11, a31);
		mcswap_var(w, a12, a32);
		mcswap_var(w, a13, a33);
		++pv;
	} else if (mc_fabs(a11) < (a21)) {
		mcswap_var(j, p0, p1);
		mcswap_var(w, a11, a21);
		mcswap_var(w, a12, a22);
		mcswap_var(w, a13, a23);
		++pv;
	}
//!# If selected pivot value != 0.
	if (a11 != 0.0) {
//!# Forming first elementary matrix.
		e121 = -(a21 / a11);
		e131 = -(a31 / a11);
	}
//!# Computing e1a=e1*a.
	e1a11 = (e111 * a11) + (e112 * a21) + (e113 * a31);
	e1a12 = (e111 * a12) + (e112 * a22) + (e113 * a32);
	e1a13 = (e111 * a13) + (e112 * a23) + (e113 * a33);

	e1a21 = (e121 * a11) + (e122 * a21) + (e123 * a31);
	e1a22 = (e121 * a12) + (e122 * a22) + (e123 * a32);
	e1a23 = (e121 * a13) + (e122 * a23) + (e123 * a33);

	e1a31 = (e131 * a11) + (e132 * a21) + (e133 * a31);
	e1a32 = (e131 * a12) + (e132 * a22) + (e133 * a32);
	e1a33 = (e131 * a13) + (e132 * a23) + (e133 * a33);

//!# Finding largest pivot bottom to top.
	if (mc_fabs(e1a32) > mc_fabs(e1a22)) {
		mcswap_var(j, p1, p2);
		mcswap_var(w, e1a21, e1a31);
		mcswap_var(w, e1a22, e1a32);
		mcswap_var(w, e1a23, e1a33);
		mcswap_var(w, e121, e131);
		++pv;
	}
//!# If selected pivot value != 0.
	if (e1a22 != 0.0) {
//!# Forming second elementary matrix e2.
		e232 = -(e1a32 / e1a22);
	}
//!# Computing U such as U=e2*e1a;
	u[0] = (e211 * e1a11) + (e212 * e1a21) + (e213 * e1a31);
	u[1] = (e211 * e1a12) + (e212 * e1a22) + (e213 * e1a32);
	u[2] = (e211 * e1a13) + (e212 * e1a23) + (e213 * e1a33);

	u[3] = (e221 * e1a11) + (e222 * e1a21) + (e223 * e1a31);
	u[4] = (e221 * e1a12) + (e222 * e1a22) + (e223 * e1a32);
	u[5] = (e221 * e1a13) + (e222 * e1a23) + (e223 * e1a33);

	u[6] = (e231 * e1a11) + (e232 * e1a21) + (e233 * e1a31);
	u[7] = (e231 * e1a12) + (e232 * e1a22) + (e233 * e1a32);
	u[8] = (e231 * e1a13) + (e232 * e1a23) + (e233 * e1a33);

//!# Computing L from e1 and e2 such as e1^-1*e2^-1.
	l[0] =  1.0;  l[1] =  0.0;  l[2] = 0.0;
	l[3] = -e121; l[4] =  1.0;  l[5] = 0.0;
	l[6] = -e131; l[7] = -e232; l[8] = 1.0;

//!# Computing P or assigning pivot indeces.
	if (wantpvi) {
		pvi[0] = p0; pvi[1] = p1; pvi[2] = p2; 
	} else if (wantp) {
		p[0] = p0 == 0 ? 1.0 : 0.0; p[1] = p0 == 1 ? 1.0 : 0.0; p[2] = p0 == 2 ? 1.0 : 0.0;
		p[3] = p1 == 0 ? 1.0 : 0.0; p[4] = p1 == 1 ? 1.0 : 0.0; p[5] = p1 == 2 ? 1.0 : 0.0;
		p[6] = p2 == 0 ? 1.0 : 0.0; p[7] = p2 == 1 ? 1.0 : 0.0; p[8] = p2 == 2 ? 1.0 : 0.0;
	}
	return pv;
}

MC_TARGET_FUNC int mc_lup3x3(const double a[9], double l[9], double u[9], double p[9], int pvi[3])
{
//!# A and L may be the same. Using a closed-form expression.
//!# Returns A=L(DU) as per Doolittle's method.
	const int wantpvi = mc_nonnullptr(pvi);
	const int wantp   = mc_nonnullptr(p);

	double w;
	int p0 = 0, p1 = 1, p2 = 2, j, pv = 0;

	double a11, a12, a13;
	double a21, a22, a23;
	double a31, a32, a33;

	double e111 =  1.0, e112 = 0.0, e113 = 0.0;
	double e121 = -0.0, e122 = 1.0, e123 = 0.0;
	double e131 = -0.0, e132 = 0.0, e133 = 1.0;

	double e211 = 1.0, e212 =  0.0, e213 = 0.0;
	double e221 = 0.0, e222 =  1.0, e223 = 0.0;
	double e231 = 0.0, e232 = -0.0, e233 = 1.0;

	double e1a11, e1a12, e1a13;
	double e1a21, e1a22, e1a23;
	double e1a31, e1a32, e1a33;

//!# Copying a.
	a11 = a[0]; a12 = a[1]; a13 = a[2];
	a21 = a[3]; a22 = a[4]; a23 = a[5];
	a31 = a[6]; a32 = a[7]; a33 = a[8];
//!# Finding largest pivot bottom to top.
	if (mc_fabs(a11) < mc_fabs(a31) && mc_fabs(a31) > mc_fabs(a21)) {
		mcswap_var(j, p0, p2);
		mcswap_var(w, a11, a31);
		mcswap_var(w, a12, a32);
		mcswap_var(w, a13, a33);
		++pv;
	} else if (mc_fabs(a11) < (a21)) {
		mcswap_var(j, p0, p1);
		mcswap_var(w, a11, a21);
		mcswap_var(w, a12, a22);
		mcswap_var(w, a13, a23);
		++pv;
	}
//!# If selected pivot value != 0.
	if (a11 != 0.0) {
//!# Forming first elementary matrix.
		e121 = -(a21 / a11);
		e131 = -(a31 / a11);
	}
//!# Computing e1a=e1*a.
	e1a11 = (e111 * a11) + (e112 * a21) + (e113 * a31);
	e1a12 = (e111 * a12) + (e112 * a22) + (e113 * a32);
	e1a13 = (e111 * a13) + (e112 * a23) + (e113 * a33);

	e1a21 = (e121 * a11) + (e122 * a21) + (e123 * a31);
	e1a22 = (e121 * a12) + (e122 * a22) + (e123 * a32);
	e1a23 = (e121 * a13) + (e122 * a23) + (e123 * a33);

	e1a31 = (e131 * a11) + (e132 * a21) + (e133 * a31);
	e1a32 = (e131 * a12) + (e132 * a22) + (e133 * a32);
	e1a33 = (e131 * a13) + (e132 * a23) + (e133 * a33);

//!# Finding largest pivot bottom to top.
	if (mc_fabs(e1a32) > mc_fabs(e1a22)) {
		mcswap_var(j, p1, p2);
		mcswap_var(w, e1a21, e1a31);
		mcswap_var(w, e1a22, e1a32);
		mcswap_var(w, e1a23, e1a33);
		mcswap_var(w, e121, e131);
		++pv;
	}
//!# If selected pivot value != 0.
	if (e1a22 != 0.0) {
//!# Forming second elementary matrix e2.
		e232 = -(e1a32 / e1a22);
	}
//!# Computing U such as U=e2*e1a;
	u[0] = (e211 * e1a11) + (e212 * e1a21) + (e213 * e1a31);
	u[1] = (e211 * e1a12) + (e212 * e1a22) + (e213 * e1a32);
	u[2] = (e211 * e1a13) + (e212 * e1a23) + (e213 * e1a33);

	u[3] = (e221 * e1a11) + (e222 * e1a21) + (e223 * e1a31);
	u[4] = (e221 * e1a12) + (e222 * e1a22) + (e223 * e1a32);
	u[5] = (e221 * e1a13) + (e222 * e1a23) + (e223 * e1a33);

	u[6] = (e231 * e1a11) + (e232 * e1a21) + (e233 * e1a31);
	u[7] = (e231 * e1a12) + (e232 * e1a22) + (e233 * e1a32);
	u[8] = (e231 * e1a13) + (e232 * e1a23) + (e233 * e1a33);

//!# Computing L from e1 and e2 such as e1^-1*e2^-1.
	l[0] =  1.0;  l[1] =  0.0;  l[2] = 0.0;
	l[3] = -e121; l[4] =  1.0;  l[5] = 0.0;
	l[6] = -e131; l[7] = -e232; l[8] = 1.0;

//!# Computing P or assigning pivot indeces.
	if (wantpvi) {
		pvi[0] = p0; pvi[1] = p1; pvi[2] = p2; 
	} else if (wantp) {
		p[0] = p0 == 0 ? 1.0 : 0.0; p[1] = p0 == 1 ? 1.0 : 0.0; p[2] = p0 == 2 ? 1.0 : 0.0;
		p[3] = p1 == 0 ? 1.0 : 0.0; p[4] = p1 == 1 ? 1.0 : 0.0; p[5] = p1 == 2 ? 1.0 : 0.0;
		p[6] = p2 == 0 ? 1.0 : 0.0; p[7] = p2 == 1 ? 1.0 : 0.0; p[8] = p2 == 2 ? 1.0 : 0.0;
	}
	return pv;
}

MC_TARGET_FUNC int mc_lup3x3l(const long double a[9], long double l[9], long double u[9], long double p[9], int pvi[3])
{
//!# A and L may be the same. Using a closed-form expression.
//!# Returns A=L(DU) as per Doolittle's method.
	const int wantpvi = mc_nonnullptr(pvi);
	const int wantp   = mc_nonnullptr(p);

	long double w;
	int p0 = 0, p1 = 1, p2 = 2, j, pv = 0;

	long double a11, a12, a13;
	long double a21, a22, a23;
	long double a31, a32, a33;

	long double e111 =  1.0L, e112 = 0.0L, e113 = 0.0L;
	long double e121 = -0.0L, e122 = 1.0L, e123 = 0.0L;
	long double e131 = -0.0L, e132 = 0.0L, e133 = 1.0L;

	long double e211 = 1.0L, e212 =  0.0L, e213 = 0.0L;
	long double e221 = 0.0L, e222 =  1.0L, e223 = 0.0L;
	long double e231 = 0.0L, e232 = -0.0L, e233 = 1.0L;

	long double e1a11, e1a12, e1a13;
	long double e1a21, e1a22, e1a23;
	long double e1a31, e1a32, e1a33;

//!# Copying a.
	a11 = a[0]; a12 = a[1]; a13 = a[2];
	a21 = a[3]; a22 = a[4]; a23 = a[5];
	a31 = a[6]; a32 = a[7]; a33 = a[8];
//!# Finding largest pivot bottom to top.
	if (mc_fabsl(a11) < mc_fabsl(a31) && mc_fabsl(a31) > mc_fabsl(a21)) {
		mcswap_var(j, p0, p2);
		mcswap_var(w, a11, a31);
		mcswap_var(w, a12, a32);
		mcswap_var(w, a13, a33);
		++pv;
	} else if (mc_fabsl(a11) < (a21)) {
		mcswap_var(j, p0, p1);
		mcswap_var(w, a11, a21);
		mcswap_var(w, a12, a22);
		mcswap_var(w, a13, a23);
		++pv;
	}
//!# If selected pivot value != 0.
	if (a11 != 0.0L) {
//!# Forming first elementary matrix.
		e121 = -(a21 / a11);
		e131 = -(a31 / a11);
	}
//!# Computing e1a=e1*a.
	e1a11 = (e111 * a11) + (e112 * a21) + (e113 * a31);
	e1a12 = (e111 * a12) + (e112 * a22) + (e113 * a32);
	e1a13 = (e111 * a13) + (e112 * a23) + (e113 * a33);

	e1a21 = (e121 * a11) + (e122 * a21) + (e123 * a31);
	e1a22 = (e121 * a12) + (e122 * a22) + (e123 * a32);
	e1a23 = (e121 * a13) + (e122 * a23) + (e123 * a33);

	e1a31 = (e131 * a11) + (e132 * a21) + (e133 * a31);
	e1a32 = (e131 * a12) + (e132 * a22) + (e133 * a32);
	e1a33 = (e131 * a13) + (e132 * a23) + (e133 * a33);

//!# Finding largest pivot bottom to top.
	if (mc_fabsl(e1a32) > mc_fabsl(e1a22)) {
		mcswap_var(j, p1, p2);
		mcswap_var(w, e1a21, e1a31);
		mcswap_var(w, e1a22, e1a32);
		mcswap_var(w, e1a23, e1a33);
		mcswap_var(w, e121, e131);
		++pv;
	}
//!# If selected pivot value != 0.
	if (e1a22 != 0.0L) {
//!# Forming second elementary matrix e2.
		e232 = -(e1a32 / e1a22);
	}
//!# Computing U such as U=e2*e1a;
	u[0] = (e211 * e1a11) + (e212 * e1a21) + (e213 * e1a31);
	u[1] = (e211 * e1a12) + (e212 * e1a22) + (e213 * e1a32);
	u[2] = (e211 * e1a13) + (e212 * e1a23) + (e213 * e1a33);

	u[3] = (e221 * e1a11) + (e222 * e1a21) + (e223 * e1a31);
	u[4] = (e221 * e1a12) + (e222 * e1a22) + (e223 * e1a32);
	u[5] = (e221 * e1a13) + (e222 * e1a23) + (e223 * e1a33);

	u[6] = (e231 * e1a11) + (e232 * e1a21) + (e233 * e1a31);
	u[7] = (e231 * e1a12) + (e232 * e1a22) + (e233 * e1a32);
	u[8] = (e231 * e1a13) + (e232 * e1a23) + (e233 * e1a33);

//!# Computing L from e1 and e2 such as e1^-1*e2^-1.
	l[0] =  1.0L; l[1] =  0.0L; l[2] = 0.0L;
	l[3] = -e121; l[4] =  1.0L; l[5] = 0.0L;
	l[6] = -e131; l[7] = -e232; l[8] = 1.0L;

//!# Computing P or assigning pivot indeces.
	if (wantpvi) {
		pvi[0] = p0; pvi[1] = p1; pvi[2] = p2; 
	} else if (wantp) {
		p[0] = p0 == 0 ? 1.0L : 0.0L; p[1] = p0 == 1 ? 1.0L : 0.0L; p[2] = p0 == 2 ? 1.0L : 0.0L;
		p[3] = p1 == 0 ? 1.0L : 0.0L; p[4] = p1 == 1 ? 1.0L : 0.0L; p[5] = p1 == 2 ? 1.0L : 0.0L;
		p[6] = p2 == 0 ? 1.0L : 0.0L; p[7] = p2 == 1 ? 1.0L : 0.0L; p[8] = p2 == 2 ? 1.0L : 0.0L;
	}
	return pv;
}

#endif /* !MC_LUP3X3_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lu3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/numa/mc_eye3x3.h>

#ifndef MC_LU3X3_H
#define MC_LU3X3_H

#pragma mark - mc_lu3x3 -

MC_TARGET_FUNC int mc_lu3x3f(const float a[9], float l[9], float u[9])
{
//!# A and L may be the same. Using a closed-form expression.
//!# Returns A=L(DU) as per Doolittle's method.
	float e1a11, e1a12, e1a13;
	float e1a21, e1a22, e1a23;
	float e1a31, e1a32, e1a33;

	float e111 =  1.0f, e112 = 0.0f, e113 = 0.0f;
	float e121 = -0.0f, e122 = 1.0f, e123 = 0.0f;
	float e131 = -0.0f, e132 = 0.0f, e133 = 1.0f;

	float e211 = 1.0f, e212 =  0.0f, e213 = 0.0f;
	float e221 = 0.0f, e222 =  1.0f, e223 = 0.0f;
	float e231 = 0.0f, e232 = -0.0f, e233 = 1.0f;

	if (a[0] != 0.0f) {
//!# Forming first elementary matrix e1.
		e121 = -(a[3] / a[0]);
		e131 = -(a[6] / a[0]);
	}
//!# Computing e1a=e1*a
	e1a11 = (e111 * a[0]) + (e112 * a[3]) + (e113 * a[6]);
	e1a12 = (e111 * a[1]) + (e112 * a[4]) + (e113 * a[7]);
	e1a13 = (e111 * a[2]) + (e112 * a[5]) + (e113 * a[8]);

	e1a21 = (e121 * a[0]) + (e122 * a[3]) + (e123 * a[6]);
	e1a22 = (e121 * a[1]) + (e122 * a[4]) + (e123 * a[7]);
	e1a23 = (e121 * a[2]) + (e122 * a[5]) + (e123 * a[8]);

	e1a31 = (e131 * a[0]) + (e132 * a[3]) + (e133 * a[6]);
	e1a32 = (e131 * a[1]) + (e132 * a[4]) + (e133 * a[7]);
	e1a33 = (e131 * a[2]) + (e132 * a[5]) + (e133 * a[8]);

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

	return 0;
}

MC_TARGET_FUNC int mc_lu3x3ff(const float a[9], double l[9], double u[9])
{
//!# Using a closed-form expression.
//!# Returns A=L(DU) as per Doolittle's method.
	double e1a11, e1a12, e1a13;
	double e1a21, e1a22, e1a23;
	double e1a31, e1a32, e1a33;

	double e111 =  1.0, e112 = 0.0, e113 = 0.0;
	double e121 = -0.0, e122 = 1.0, e123 = 0.0;
	double e131 = -0.0, e132 = 0.0, e133 = 1.0;

	double e211 = 1.0, e212 =  0.0, e213 = 0.0;
	double e221 = 0.0, e222 =  1.0, e223 = 0.0;
	double e231 = 0.0, e232 = -0.0, e233 = 1.0;

	if (mc_cast(double, a[0]) != 0.0) {
//!# Forming first elementary matrix e1.
		e121 = -(mc_cast(double, a[3]) / mc_cast(const double, a[0]));
		e131 = -(mc_cast(double, a[6]) / mc_cast(const double, a[0]));
	}
//!# Computing e1a=e1*a
	e1a11 = (e111 * mc_cast(double, a[0])) + (e112 * mc_cast(double, a[3])) + (e113 * mc_cast(double, a[6]));
	e1a12 = (e111 * mc_cast(double, a[1])) + (e112 * mc_cast(double, a[4])) + (e113 * mc_cast(double, a[7]));
	e1a13 = (e111 * mc_cast(double, a[2])) + (e112 * mc_cast(double, a[5])) + (e113 * mc_cast(double, a[8]));

	e1a21 = (e121 * mc_cast(double, a[0])) + (e122 * mc_cast(double, a[3])) + (e123 * mc_cast(double, a[6]));
	e1a22 = (e121 * mc_cast(double, a[1])) + (e122 * mc_cast(double, a[4])) + (e123 * mc_cast(double, a[7]));
	e1a23 = (e121 * mc_cast(double, a[2])) + (e122 * mc_cast(double, a[5])) + (e123 * mc_cast(double, a[8]));

	e1a31 = (e131 * mc_cast(double, a[0])) + (e132 * mc_cast(double, a[3])) + (e133 * mc_cast(double, a[6]));
	e1a32 = (e131 * mc_cast(double, a[1])) + (e132 * mc_cast(double, a[4])) + (e133 * mc_cast(double, a[7]));
	e1a33 = (e131 * mc_cast(double, a[2])) + (e132 * mc_cast(double, a[5])) + (e133 * mc_cast(double, a[8]));

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

	return 0;
}

MC_TARGET_FUNC int mc_lu3x3(const double a[9], double l[9], double u[9])
{
//!# A and L may be the same. Using a closed-form expression.
//!# Returns A=L(DU) as per Doolittle's method.
	double e1a11, e1a12, e1a13;
	double e1a21, e1a22, e1a23;
	double e1a31, e1a32, e1a33;

	double e111 =  1.0, e112 = 0.0, e113 = 0.0;
	double e121 = -0.0, e122 = 1.0, e123 = 0.0;
	double e131 = -0.0, e132 = 0.0, e133 = 1.0;

	double e211 = 1.0, e212 =  0.0, e213 = 0.0;
	double e221 = 0.0, e222 =  1.0, e223 = 0.0;
	double e231 = 0.0, e232 = -0.0, e233 = 1.0;

	if (a[0] != 0.0) {
//!# Forming first elementary matrix e1.
		e121 = -(a[3] / a[0]);
		e131 = -(a[6] / a[0]);
	}
//!# Computing e1a=e1*a
	e1a11 = (e111 * a[0]) + (e112 * a[3]) + (e113 * a[6]);
	e1a12 = (e111 * a[1]) + (e112 * a[4]) + (e113 * a[7]);
	e1a13 = (e111 * a[2]) + (e112 * a[5]) + (e113 * a[8]);

	e1a21 = (e121 * a[0]) + (e122 * a[3]) + (e123 * a[6]);
	e1a22 = (e121 * a[1]) + (e122 * a[4]) + (e123 * a[7]);
	e1a23 = (e121 * a[2]) + (e122 * a[5]) + (e123 * a[8]);

	e1a31 = (e131 * a[0]) + (e132 * a[3]) + (e133 * a[6]);
	e1a32 = (e131 * a[1]) + (e132 * a[4]) + (e133 * a[7]);
	e1a33 = (e131 * a[2]) + (e132 * a[5]) + (e133 * a[8]);

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

	return 0;
}

MC_TARGET_FUNC int mc_lu3x3l(const long double a[9], long double l[9], long double u[9])
{
//!# A and L may be the same. Using a closed-form expression.
//!# Returns A=L(DU) as per Doolittle's method.
	long double e1a11, e1a12, e1a13;
	long double e1a21, e1a22, e1a23;
	long double e1a31, e1a32, e1a33;

	long double e111 =  1.0L, e112 = 0.0L, e113 = 0.0L;
	long double e121 = -0.0L, e122 = 1.0L, e123 = 0.0L;
	long double e131 = -0.0L, e132 = 0.0L, e133 = 1.0L;

	long double e211 = 1.0L, e212 =  0.0L, e213 = 0.0L;
	long double e221 = 0.0L, e222 =  1.0L, e223 = 0.0L;
	long double e231 = 0.0L, e232 = -0.0L, e233 = 1.0L;

	if (a[0] != 0.0L) {
//!# Forming first elementary matrix e1.
		e121 = -(a[3] / a[0]);
		e131 = -(a[6] / a[0]);
	}
//!# Computing e1a=e1*a
	e1a11 = (e111 * a[0]) + (e112 * a[3]) + (e113 * a[6]);
	e1a12 = (e111 * a[1]) + (e112 * a[4]) + (e113 * a[7]);
	e1a13 = (e111 * a[2]) + (e112 * a[5]) + (e113 * a[8]);

	e1a21 = (e121 * a[0]) + (e122 * a[3]) + (e123 * a[6]);
	e1a22 = (e121 * a[1]) + (e122 * a[4]) + (e123 * a[7]);
	e1a23 = (e121 * a[2]) + (e122 * a[5]) + (e123 * a[8]);

	e1a31 = (e131 * a[0]) + (e132 * a[3]) + (e133 * a[6]);
	e1a32 = (e131 * a[1]) + (e132 * a[4]) + (e133 * a[7]);
	e1a33 = (e131 * a[2]) + (e132 * a[5]) + (e133 * a[8]);

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

	return 0;
}

#endif /* !MC_LU3X3_H */

/* EOF */
//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_svd3x3.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/numa/mc_eigsy3x3.h>
#include <macadam/details/numa/mc_eye3x3.h>
#include <macadam/details/numa/mc_mulab3x3.h>
#include <macadam/details/numa/mc_mulatb3x3.h>
#include <macadam/details/numa/mc_qrgv3x3.h>
#include <macadam/mcswap.h>

#ifndef MC_SVD3X3_H
#define MC_SVD3X3_H

#pragma mark - mc_svd3x3 -

MC_TARGET_FUNC int mc_svd3x3f(const float a[9], float u[9], float s[9], float v[9])
{
//!# The main result SVD provides is that we can write an m by n matrix A
//!# such as U'*A=S*V' with:
//!#     - U is an [m x p] orthogonal matrix. The left-singular vectors of A are a set of orthonormal eigenvectors of AA'.
//!#     - S is an [n x p] diagonal matrix. The non-negative singular values of A (found on the diagonal entries of S) are
//!#       the square roots of the non-negative eigenvalues of both AA' and A'A.
//!#     - V is an [p x p] orthogonal matrix. The right-singular vectors of A are a set of orthonormal eigenvectors of A'A.
//!#     - p=min(m, n) and in this particular case we have m=3, n=3 hence p=3.
	float w;

//!# Step 1: Forming A'*A storing temporarily result into U.
	mc_mulatb3x3f(u, a, a);

//!# Step 2: Computing V of A'*A i.e right-singular vectors.
	if (0 == mc_eigsy3x3f(u, s, v)) {
		mc_eye3x3f(s, 0);
		mc_mulab3x3f(u, a, v);
//!# Step 3: Computing singular-values and U such as U=A*V*S^-1.
		if (0 == mc_qrgv3x3f(u, u, s)) {
//!# Step 4: Decimeting off-diagonal residual elements.
			s[1] = 0.0f; s[2] = 0.0f;
			s[3] = 0.0f; s[5] = 0.0f;
			s[6] = 0.0f; s[7] = 0.0f;
//!# Step 5: Changing sign if required.
			if (mc_copysignf(1.0f, s[0]) < 0.0f) {
				s[0] = -s[0];
				v[0] = -v[0];
				v[3] = -v[3];
				v[6] = -v[6];
			}
			if (mc_copysignf(1.0f, s[4]) < 0.0f) {
				s[4] = -s[4];
				v[1] = -v[1];
				v[4] = -v[4];
				v[7] = -v[7];
			}
			if (mc_copysignf(1.0f, s[8]) < 0.0f) {
				s[8] = -s[8];
				v[2] = -v[2];
				v[5] = -v[5];
				v[8] = -v[8];
			}
//!# Step 6: Reordering singular-values and basis (descending i.e largest first).
			if (s[0] < s[4]) {
				mcswap_var(w, s[0], s[4]);

				mcswap_var(w, v[0], v[1]);
				mcswap_var(w, v[3], v[4]);
				mcswap_var(w, v[6], v[7]);

				mcswap_var(w, u[0], u[1]);
				mcswap_var(w, u[3], u[4]);
				mcswap_var(w, u[6], u[7]);
			}
			if (s[0] < s[8]) {
				mcswap_var(w, s[0], s[8]);

				mcswap_var(w, v[0], v[2]);
				mcswap_var(w, v[3], v[5]);
				mcswap_var(w, v[6], v[8]);

				mcswap_var(w, u[0], u[2]);
				mcswap_var(w, u[3], u[5]);
				mcswap_var(w, u[6], u[8]);
			}
			if (s[4] < s[8]) {
				mcswap_var(w, s[4], s[8]);

				mcswap_var(w, v[1], v[2]);
				mcswap_var(w, v[4], v[5]);
				mcswap_var(w, v[7], v[8]);

				mcswap_var(w, u[1], u[2]);
				mcswap_var(w, u[4], u[5]);
				mcswap_var(w, u[7], u[8]);
			}
			return 0;
		}
		return -2;
	}
	return -1;
}

MC_TARGET_FUNC int mc_svd3x3ff(const float a[9], double u[9], double s[9], double v[9])
{
//!# The main result SVD provides is that we can write an m by n matrix A
//!# such as U'*A=S*V' with:
//!#     - U is an [m x p] orthogonal matrix. The left-singular vectors of A are a set of orthonormal eigenvectors of AA'.
//!#     - S is an [n x p] diagonal matrix. The non-negative singular values of A (found on the diagonal entries of S) are
//!#       the square roots of the non-negative eigenvalues of both AA' and A'A.
//!#     - V is an [p x p] orthogonal matrix. The right-singular vectors of A are a set of orthonormal eigenvectors of A'A.
//!#     - p=min(m, n) and in this particular case we have m=3, n=3 hence p=3.
	double w;

//!# Step 1: Forming A'*A storing temporarily result into U.
	mc_mulatb3x3ff(u, a, a);

//!# Step 2: Computing V of A'*A i.e right-singular vectors.
	if (0 == mc_eigsy3x3(u, s, v)) {
		mc_eye3x3(s, 0);
		mc_mulab3x3fd(u, a, v);
//!# Step 3: Computing singular-values and U such as U=A*V*S^-1.
		if (0 == mc_qrgv3x3(u, u, s)) {
//!# Step 4: Decimeting off-diagonal residual elements.
			s[1] = 0.0; s[2] = 0.0;
			s[3] = 0.0; s[5] = 0.0;
			s[6] = 0.0; s[7] = 0.0;
//!# Step 5: Changing sign if required.
			if (mc_copysign(1.0, s[0]) < 0.0) {
				s[0] = -s[0];
				v[0] = -v[0];
				v[3] = -v[3];
				v[6] = -v[6];
			}
			if (mc_copysign(1.0, s[4]) < 0.0) {
				s[4] = -s[4];
				v[1] = -v[1];
				v[4] = -v[4];
				v[7] = -v[7];
			}
			if (mc_copysign(1.0, s[8]) < 0.0) {
				s[8] = -s[8];
				v[2] = -v[2];
				v[5] = -v[5];
				v[8] = -v[8];
			}
//!# Step 6: Reordering singular-values and basis (descending i.e largest first).
			if (s[0] < s[4]) {
				mcswap_var(w, s[0], s[4]);

				mcswap_var(w, v[0], v[1]);
				mcswap_var(w, v[3], v[4]);
				mcswap_var(w, v[6], v[7]);

				mcswap_var(w, u[0], u[1]);
				mcswap_var(w, u[3], u[4]);
				mcswap_var(w, u[6], u[7]);
			}
			if (s[0] < s[8]) {
				mcswap_var(w, s[0], s[8]);

				mcswap_var(w, v[0], v[2]);
				mcswap_var(w, v[3], v[5]);
				mcswap_var(w, v[6], v[8]);

				mcswap_var(w, u[0], u[2]);
				mcswap_var(w, u[3], u[5]);
				mcswap_var(w, u[6], u[8]);
			}
			if (s[4] < s[8]) {
				mcswap_var(w, s[4], s[8]);

				mcswap_var(w, v[1], v[2]);
				mcswap_var(w, v[4], v[5]);
				mcswap_var(w, v[7], v[8]);

				mcswap_var(w, u[1], u[2]);
				mcswap_var(w, u[4], u[5]);
				mcswap_var(w, u[7], u[8]);
			}
			return 0;
		}
		return -2;
	}
	return -1;
}

MC_TARGET_FUNC int mc_svd3x3(const double a[9], double u[9], double s[9], double v[9])
{
//!# The main result SVD provides is that we can write an m by n matrix A
//!# such as U'*A=S*V' with:
//!#     - U is an [m x p] orthogonal matrix. The left-singular vectors of A are a set of orthonormal eigenvectors of AA'.
//!#     - S is an [n x p] diagonal matrix. The non-negative singular values of A (found on the diagonal entries of S) are
//!#       the square roots of the non-negative eigenvalues of both AA' and A'A.
//!#     - V is an [p x p] orthogonal matrix. The right-singular vectors of A are a set of orthonormal eigenvectors of A'A.
//!#     - p=min(m, n) and in this particular case we have m=3, n=3 hence p=3.
	double w;

//!# Step 1: Forming A'*A storing temporarily result into U.
	mc_mulatb3x3(u, a, a);

//!# Step 2: Computing V of A'*A i.e right-singular vectors.
	if (0 == mc_eigsy3x3(u, s, v)) {
		mc_eye3x3(s, 0);
		mc_mulab3x3(u, a, v);
//!# Step 3: Computing singular-values and U such as U=A*V*S^-1.
		if (0 == mc_qrgv3x3(u, u, s)) {
//!# Step 4: Decimeting off-diagonal residual elements.
			s[1] = 0.0; s[2] = 0.0;
			s[3] = 0.0; s[5] = 0.0;
			s[6] = 0.0; s[7] = 0.0;
//!# Step 5: Changing sign if required.
			if (mc_copysign(1.0, s[0]) < 0.0) {
				s[0] = -s[0];
				v[0] = -v[0];
				v[3] = -v[3];
				v[6] = -v[6];
			}
			if (mc_copysign(1.0, s[4]) < 0.0) {
				s[4] = -s[4];
				v[1] = -v[1];
				v[4] = -v[4];
				v[7] = -v[7];
			}
			if (mc_copysign(1.0, s[8]) < 0.0) {
				s[8] = -s[8];
				v[2] = -v[2];
				v[5] = -v[5];
				v[8] = -v[8];
			}
//!# Step 6: Reordering singular-values and basis (descending i.e largest first).
			if (s[0] < s[4]) {
				mcswap_var(w, s[0], s[4]);

				mcswap_var(w, v[0], v[1]);
				mcswap_var(w, v[3], v[4]);
				mcswap_var(w, v[6], v[7]);

				mcswap_var(w, u[0], u[1]);
				mcswap_var(w, u[3], u[4]);
				mcswap_var(w, u[6], u[7]);
			}
			if (s[0] < s[8]) {
				mcswap_var(w, s[0], s[8]);

				mcswap_var(w, v[0], v[2]);
				mcswap_var(w, v[3], v[5]);
				mcswap_var(w, v[6], v[8]);

				mcswap_var(w, u[0], u[2]);
				mcswap_var(w, u[3], u[5]);
				mcswap_var(w, u[6], u[8]);
			}
			if (s[4] < s[8]) {
				mcswap_var(w, s[4], s[8]);

				mcswap_var(w, v[1], v[2]);
				mcswap_var(w, v[4], v[5]);
				mcswap_var(w, v[7], v[8]);

				mcswap_var(w, u[1], u[2]);
				mcswap_var(w, u[4], u[5]);
				mcswap_var(w, u[7], u[8]);
			}
			return 0;
		}
		return -2;
	}
	return -1;
}

MC_TARGET_FUNC int mc_svd3x3l(const long double a[9], long double u[9], long double s[9], long double v[9])
{
//!# The main result SVD provides is that we can write an m by n matrix A
//!# such as U'*A=S*V' with:
//!#     - U is an [m x p] orthogonal matrix. The left-singular vectors of A are a set of orthonormal eigenvectors of AA'.
//!#     - S is an [n x p] diagonal matrix. The non-negative singular values of A (found on the diagonal entries of S) are
//!#       the square roots of the non-negative eigenvalues of both AA' and A'A.
//!#     - V is an [p x p] orthogonal matrix. The right-singular vectors of A are a set of orthonormal eigenvectors of A'A.
//!#     - p=min(m, n) and in this particular case we have m=3, n=3 hence p=3.
	long double w;

//!# Step 1: Forming A'*A storing temporarily result into U.
	mc_mulatb3x3l(u, a, a);

//!# Step 2: Computing V of A'*A i.e right-singular vectors.
	if (0 == mc_eigsy3x3l(u, s, v)) {
		mc_eye3x3l(s, 0);
		mc_mulab3x3l(u, a, v);
//!# Step 3: Computing singular-values and U such as U=A*V*S^-1.
		if (0 == mc_qrgv3x3l(u, u, s)) {
//!# Step 4: Decimeting off-diagonal residual elements.
			s[1] = 0.0L; s[2] = 0.0L;
			s[3] = 0.0L; s[5] = 0.0L;
			s[6] = 0.0L; s[7] = 0.0L;
//!# Step 5: Changing sign if required.
			if (mc_copysignl(1.0L, s[0]) < 0.0L) {
				s[0] = -s[0];
				v[0] = -v[0];
				v[3] = -v[3];
				v[6] = -v[6];
			}
			if (mc_copysignl(1.0L, s[4]) < 0.0L) {
				s[4] = -s[4];
				v[1] = -v[1];
				v[4] = -v[4];
				v[7] = -v[7];
			}
			if (mc_copysignl(1.0L, s[8]) < 0.0L) {
				s[8] = -s[8];
				v[2] = -v[2];
				v[5] = -v[5];
				v[8] = -v[8];
			}
//!# Step 6: Reordering singular-values and basis (descending i.e largest first).
			if (s[0] < s[4]) {
				mcswap_var(w, s[0], s[4]);

				mcswap_var(w, v[0], v[1]);
				mcswap_var(w, v[3], v[4]);
				mcswap_var(w, v[6], v[7]);

				mcswap_var(w, u[0], u[1]);
				mcswap_var(w, u[3], u[4]);
				mcswap_var(w, u[6], u[7]);
			}
			if (s[0] < s[8]) {
				mcswap_var(w, s[0], s[8]);

				mcswap_var(w, v[0], v[2]);
				mcswap_var(w, v[3], v[5]);
				mcswap_var(w, v[6], v[8]);

				mcswap_var(w, u[0], u[2]);
				mcswap_var(w, u[3], u[5]);
				mcswap_var(w, u[6], u[8]);
			}
			if (s[4] < s[8]) {
				mcswap_var(w, s[4], s[8]);

				mcswap_var(w, v[1], v[2]);
				mcswap_var(w, v[4], v[5]);
				mcswap_var(w, v[7], v[8]);

				mcswap_var(w, u[1], u[2]);
				mcswap_var(w, u[4], u[5]);
				mcswap_var(w, u[7], u[8]);
			}
			return 0;
		}
		return -2;
	}
	return -1;
}

#endif /* !MC_SVD3X3_H */

/* EOF */
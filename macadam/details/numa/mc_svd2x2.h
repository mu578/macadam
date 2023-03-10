//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_svd2x2.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_copysign.h>
#include <macadam/details/numa/mc_eigsy2x2.h>
#include <macadam/details/numa/mc_eye2x2.h>
#include <macadam/details/numa/mc_mulab2x2.h>
#include <macadam/details/numa/mc_mulatb2x2.h>
#include <macadam/details/numa/mc_qrgs2x2.h>
#include <macadam/mcswap.h>

#ifndef MC_SVD2X2_H
#define MC_SVD2X2_H

#pragma mark - mc_svd2x2 -

MC_TARGET_FUNC int mc_svd2x2f(const float a[4], float u[4], float s[4], float v[4])
{
//!# The main result SVD provides is that we can write an m by n matrix A
//!# such as U'*A=S*V' with:
//!#     - U is an [m x p] orthogonal matrix. The left-singular vectors of A are a set of orthonormal eigenvectors of AA'.
//!#     - S is an [n x p] diagonal matrix. The non-negative singular values of A (found on the diagonal entries of S) are
//!#       the square roots of the non-negative eigenvalues of both AA' and A'A.
//!#     - V is an [p x p] orthogonal matrix. The right-singular vectors of A are a set of orthonormal eigenvectors of A'A.
//!#     - p=min(m, n) and in this particular case we have m=2, n=2 hence p=2.
	float w;

//!# Step 1: Forming A'*A storing temporarily result into U.
	mc_mulatb2x2f(u, a, a);

//!# Step 2: Computing V of A'*A i.e right-singular vectors.
	if (0 == mc_eigsy2x2f(u, s, v)) {
		mc_eye2x2f(s, 0);
		mc_mulab2x2f(u, a, v);
//!# Step 3: Computing singular-values and U such as U=A*V*S^-1.
		if (0 == mc_qrgs2x2f(u, u, s)) {
//!# Step 4: Decimeting off-diagonal residual elements.
			s[1] = 0.0f; s[2] = 0.0f;
//!# Step 5: Changing sign if required.
			if (mc_copysignf(1.0f, s[0]) < 0.0f) {
				s[0] = -s[0];

				v[0] = -v[0];
				v[2] = -v[2];
			}
			if (mc_copysignf(1.0f, s[3]) < 0.0f) {
				s[3] = -s[3];

				v[1] = -v[1];
				v[3] = -v[3];
			}
//!# Step 6: Reordering singular-values and basis (descending i.e largest first).
			if (s[0] < s[3]) {
				mcswap_var(w, s[0], s[3]);

				mcswap_var(w, v[0], v[1]);
				mcswap_var(w, v[2], v[3]);

				mcswap_var(w, u[0], u[1]);
				mcswap_var(w, u[2], u[3]);
			}
			return 0;
		}
		return -2;
	}
	return -1;
}

MC_TARGET_FUNC int mc_svd2x2ff(const float a[4], double u[4], double s[4], double v[4])
{
//!# The main result SVD provides is that we can write an m by n matrix A
//!# such as U'*A=S*V' with:
//!#     - U is an [m x p] orthogonal matrix. The left-singular vectors of A are a set of orthonormal eigenvectors of AA'.
//!#     - S is an [n x p] diagonal matrix. The non-negative singular values of A (found on the diagonal entries of S) are
//!#       the square roots of the non-negative eigenvalues of both AA' and A'A.
//!#     - V is an [p x p] orthogonal matrix. The right-singular vectors of A are a set of orthonormal eigenvectors of A'A.
//!#     - p=min(m, n) and in this particular case we have m=2, n=2 hence p=2.
	double w;

//!# Step 1: Forming A'*A storing temporarily result into U.
	mc_mulatb2x2ff(u, a, a);

//!# Step 2: Computing V of A'*A i.e right-singular vectors.
	if (0 == mc_eigsy2x2(u, s, v)) {
		mc_eye2x2(s, 0);
		mc_mulab2x2fd(u, a, v);
//!# Step 3: Computing singular-values and U such as U=A*V*S^-1.
		if (0 == mc_qrgs2x2(u, u, s)) {
//!# Step 4: Decimeting off-diagonal residual elements.
			s[1] = 0.0; s[2] = 0.0;
//!# Step 5: Changing sign if required.
			if (mc_copysign(1.0, s[0]) < 0.0) {
				s[0] = -s[0];

				v[0] = -v[0];
				v[2] = -v[2];
			}
			if (mc_copysign(1.0, s[3]) < 0.0) {
				s[3] = -s[3];

				v[1] = -v[1];
				v[3] = -v[3];
			}
//!# Step 6: Reordering singular-values and basis (descending i.e largest first).
			if (s[0] < s[3]) {
				mcswap_var(w, s[0], s[3]);

				mcswap_var(w, v[0], v[1]);
				mcswap_var(w, v[2], v[3]);

				mcswap_var(w, u[0], u[1]);
				mcswap_var(w, u[2], u[3]);
			}
			return 0;
		}
		return -2;
	}
	return -1;
}

MC_TARGET_FUNC int mc_svd2x2(const double a[4], double u[4], double s[4], double v[4])
{
//!# The main result SVD provides is that we can write an m by n matrix A
//!# such as U'*A=S*V' with:
//!#     - U is an [m x p] orthogonal matrix. The left-singular vectors of A are a set of orthonormal eigenvectors of AA'.
//!#     - S is an [n x p] diagonal matrix. The non-negative singular values of A (found on the diagonal entries of S) are
//!#       the square roots of the non-negative eigenvalues of both AA' and A'A.
//!#     - V is an [p x p] orthogonal matrix. The right-singular vectors of A are a set of orthonormal eigenvectors of A'A.
//!#     - p=min(m, n) and in this particular case we have m=2, n=2 hence p=2.
	double w;

//!# Step 1: Forming A'*A storing temporarily result into U.
	mc_mulatb2x2(u, a, a);

//!# Step 2: Computing V of A'*A i.e right-singular vectors.
	if (0 == mc_eigsy2x2(u, s, v)) {
		mc_eye2x2(s, 0);
		mc_mulab2x2(u, a, v);
//!# Step 3: Computing singular-values and U such as U=A*V*S^-1.
		if (0 == mc_qrgs2x2(u, u, s)) {
//!# Step 4: Decimeting off-diagonal residual elements.
			s[1] = 0.0; s[2] = 0.0;
//!# Step 5: Changing sign if required.
			if (mc_copysign(1.0, s[0]) < 0.0) {
				s[0] = -s[0];

				v[0] = -v[0];
				v[2] = -v[2];
			}
			if (mc_copysign(1.0, s[3]) < 0.0) {
				s[3] = -s[3];

				v[1] = -v[1];
				v[3] = -v[3];
			}
//!# Step 6: Reordering singular-values and basis (descending i.e largest first).
			if (s[0] < s[3]) {
				mcswap_var(w, s[0], s[3]);

				mcswap_var(w, v[0], v[1]);
				mcswap_var(w, v[2], v[3]);

				mcswap_var(w, u[0], u[1]);
				mcswap_var(w, u[2], u[3]);
			}
			return 0;
		}
		return -2;
	}
	return -1;
}

MC_TARGET_FUNC int mc_svd2x2l(const long double a[4], long double u[4], long double s[4], long double v[4])
{
//!# The main result SVD provides is that we can write an m by n matrix A
//!# such as U'*A=S*V' with:
//!#     - U is an [m x p] orthogonal matrix. The left-singular vectors of A are a set of orthonormal eigenvectors of AA'.
//!#     - S is an [n x p] diagonal matrix. The non-negative singular values of A (found on the diagonal entries of S) are
//!#       the square roots of the non-negative eigenvalues of both AA' and A'A.
//!#     - V is an [p x p] orthogonal matrix. The right-singular vectors of A are a set of orthonormal eigenvectors of A'A.
//!#     - p=min(m, n) and in this particular case we have m=2, n=2 hence p=2.
	long double w;

//!# Step 1: Forming A'*A storing temporarily result into U.
	mc_mulatb2x2l(u, a, a);

//!# Step 2: Computing V of A'*A i.e right-singular vectors.
	if (0 == mc_eigsy2x2l(u, s, v)) {
		mc_eye2x2l(s, 0);
		mc_mulab2x2l(u, a, v);
//!# Step 3: Computing singular-values and U such as U=A*V*S^-1.
		if (0 == mc_qrgs2x2l(u, u, s)) {
//!# Step 4: Decimeting off-diagonal residual elements.
			s[1] = 0.0L; s[2] = 0.0L;
//!# Step 5: Changing sign if required.
			if (mc_copysignl(1.0L, s[0]) < 0.0L) {
				s[0] = -s[0];

				v[0] = -v[0];
				v[2] = -v[2];
			}
			if (mc_copysignl(1.0L, s[3]) < 0.0L) {
				s[3] = -s[3];

				v[1] = -v[1];
				v[3] = -v[3];
			}
//!# Step 6: Reordering singular-values and basis (descending i.e largest first).
			if (s[0] < s[3]) {
				mcswap_var(w, s[0], s[3]);

				mcswap_var(w, v[0], v[1]);
				mcswap_var(w, v[2], v[3]);

				mcswap_var(w, u[0], u[1]);
				mcswap_var(w, u[2], u[3]);
			}
			return 0;
		}
		return -2;
	}
	return -1;
}

#endif /* !MC_SVD2X2_H */

/* EOF */
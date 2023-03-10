//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_legendre_qn.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/math/mc_acoth.h>
#include <macadam/details/math/mc_atanh.h>
#include <macadam/details/math/mc_fabs.h>
#include <macadam/details/math/mc_isinf.h>
#include <macadam/details/math/mc_isnan.h>
#include <macadam/details/math/mc_legendre_pn.h>
#include <macadam/details/math/mc_log1p.h>
#include <macadam/details/math/mc_raise2.h>
#include <macadam/details/math/mc_raise3.h>
#include <macadam/details/math/mc_raise4.h>
#include <macadam/details/math/mc_rsqr.h>

#ifndef MC_LEGENDRE_QN_H
#define MC_LEGENDRE_QN_H

#pragma mark - mc_legendre_q0 -

MC_TARGET_PROC float mc_legendre_q0f(const float x)
{
//!# Legendre functions of the second kind, degree 0.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	if (mc_fabsf(x) < 1.0f) {
		return mc_atanhf(x);
	}
	return mc_acothf(x);
}

MC_TARGET_PROC double mc_legendre_q0(const double x)
{
//!# Legendre functions of the second kind, degree 0.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	if (mc_fabs(x) < 1.0) {
		return mc_atanh(x);
	}
	return mc_acoth(x);
}

MC_TARGET_PROC long double mc_legendre_q0l(const long double x)
{
//!# Legendre functions of the second kind, degree 0.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	if (mc_fabsl(x) < 1.0L) {
		return mc_atanhl(x);
	}
	return mc_acothl(x);
}

#pragma mark - mc_legendre_q1 -

MC_TARGET_PROC float mc_legendre_q1f(const float x)
{
//!# Legendre functions of the second kind, degree 1.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	float q1 = 0.0f;
	float w, y, f, s, k;

	if (mc_isnan(x) || mc_isinf(x)) {
		return MCK_NAN;
	}
	w = mc_fabsf(x);
	if (w < 1.0f) {
		q1 = mc_atanhf(x) * x - 1.0f;
	} else if (w < 2.0f) {
		q1 = 0.5f * mc_log1pf(2.0f / (w - 1.0f)) * w - 1.0f;
	} else {
		y = mc_rsqrf(x); f = 3.0f; s = 1.0f; k = 5.0f;
		do {
			f = f * y;
			w = f / k;
			s = s + w;
			k = k + 2.0f;
		} while (w <= s * MCLIMITS_EPSILONF);
		q1 = (s * y) * MCK_KF(MCK_1_3);
	}
	return q1;
}

MC_TARGET_PROC double mc_legendre_q1(const double x)
{
//!# Legendre functions of the second kind, degree 1.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	double q1 = 0.0;
	double w, y, f, s, k;

	if (mc_isnan(x) || mc_isinf(x)) {
		return MCK_NAN;
	}
	w = mc_fabs(x);
	if (w < 1.0) {
		q1 = mc_atanh(x) * x - 1.0;
	} else if (w < 2.0) {
		q1 = 0.5f * mc_log1p(2.0 / (w - 1.0)) * w - 1.0;
	} else {
		y = mc_rsqr(x); f = 3.0; s = 1.0; k = 5.0;
		do {
			f = f * y;
			w = f / k;
			s = s + w;
			k = k + 2.0;
		} while (w <= s * MCLIMITS_EPSILON);
		q1 = (s * y) * MCK_K(MCK_1_3);
	}
	return q1;
}

MC_TARGET_PROC long double mc_legendre_q1l(const long double x)
{
//!# Legendre functions of the second kind, degree 1.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	long double q1 = 0.0L;
	long double w, y, f, s, k;

	if (mc_isnan(x) || mc_isinf(x)) {
		return MCK_NAN;
	}
	w = mc_fabsl(x);
	if (w < 1.0L) {
		q1 = mc_atanhl(x) * x - 1.0L;
	} else if (w < 2.0L) {
		q1 = 0.5f * mc_log1pl(2.0L / (w - 1.0L)) * w - 1.0L;
	} else {
		y = mc_rsqrl(x); f = 3.0L; s = 1.0L; k = 5.0L;
		do {
			f = f * y;
			w = f / k;
			s = s + w;
			k = k + 2.0L;
		} while (w <= s * MCLIMITS_EPSILONL);
		q1 = (s * y) * MCK_KL(MCK_1_3);
	}
	return q1;
}

#pragma mark - mc_legendre_q2 -

MC_TARGET_PROC float mc_legendre_q2f(const float x)
{
//!# Legendre functions of the second kind, degree 2.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	const float p2 = mc_legendre_p2f(x);
	const float q0 = mc_legendre_q0f(x);
	return p2 * q0 - (MCK_KF(MCK_3_2) * x);
}

MC_TARGET_PROC double mc_legendre_q2(const double x)
{
//!# Legendre functions of the second kind, degree 2.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	const double p2 = mc_legendre_p2(x);
	const double q0 = mc_legendre_q0(x);
	return p2 * q0 - (MCK_K(MCK_3_2) * x);
}

MC_TARGET_PROC long double mc_legendre_q2l(const long double x)
{
//!# Legendre functions of the second kind, degree 2.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	const long double p2 = mc_legendre_p2l(x);
	const long double q0 = mc_legendre_q0l(x);
	return p2 * q0 - (MCK_KL(MCK_3_2) * x);
}

#pragma mark - mc_legendre_q3 -

MC_TARGET_PROC float mc_legendre_q3f(const float x)
{
//!# Legendre functions of the second kind, degree 3.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	const float p3 = mc_legendre_p3f(x);
	const float q0 = mc_legendre_q0f(x);
	return p3 * q0 - (MCK_KF(MCK_5_2) * mc_raise2f(x)) + MCK_KF(MCK_2_3);
}

MC_TARGET_PROC double mc_legendre_q3(const double x)
{
//!# Legendre functions of the second kind, degree 3.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	const double p3 = mc_legendre_p3(x);
	const double q0 = mc_legendre_q0(x);
	return p3 * q0 - (MCK_K(MCK_5_2) * mc_raise2(x)) + MCK_K(MCK_2_3);
}

MC_TARGET_PROC long double mc_legendre_q3l(const long double x)
{
//!# Legendre functions of the second kind, degree 3.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	const long double p3 = mc_legendre_p3l(x);
	const long double q0 = mc_legendre_q0l(x);
	return p3 * q0 - (MCK_KL(MCK_5_2) * mc_raise2l(x)) + MCK_KL(MCK_2_3);
}

#pragma mark - mc_legendre_q4 -

MC_TARGET_PROC float mc_legendre_q4f(const float x)
{
//!# Legendre functions of the second kind, degree 4.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	const float p4 = mc_legendre_p4f(x);
	const float q0 = mc_legendre_q0f(x);
	return p4 * q0 - (MCK_KF(MCK_35_8) * mc_raise3f(x)) + (MCK_KF(MCK_55_24) * x);
}

MC_TARGET_PROC double mc_legendre_q4(const double x)
{
//!# Legendre functions of the second kind, degree 4.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	const double p4 = mc_legendre_p4(x);
	const double q0 = mc_legendre_q0(x);
	return p4 * q0 - (MCK_K(MCK_35_8) * mc_raise3(x)) + (MCK_K(MCK_55_24) * x);
}

MC_TARGET_PROC long double mc_legendre_q4l(const long double x)
{
//!# Legendre functions of the second kind, degree 4.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	const long double p4 = mc_legendre_p4l(x);
	const long double q0 = mc_legendre_q0l(x);
	return p4 * q0 - (MCK_KL(MCK_35_8) * mc_raise3l(x)) + (MCK_KL(MCK_55_24) * x);
}

#pragma mark - mc_legendre_q5 -

MC_TARGET_PROC float mc_legendre_q5f(const float x)
{
//!# Legendre functions of the second kind, degree 5.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	const float p5 = mc_legendre_p5f(x);
	const float q0 = mc_legendre_q0f(x);
	return p5 * q0 - (MCK_KF(MCK_63_8) * mc_raise4f(x)) + (MCK_KF(MCK_49_8) * mc_raise2f(x)) - MCK_KF(MCK_8_15);
}

MC_TARGET_PROC double mc_legendre_q5(const double x)
{
//!# Legendre functions of the second kind, degree 5.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	const double p5 = mc_legendre_p5(x);
	const double q0 = mc_legendre_q0(x);
	return p5 * q0 - (MCK_K(MCK_63_8) * mc_raise4(x)) + (MCK_K(MCK_49_8) * mc_raise2(x)) - MCK_K(MCK_8_15);
}

MC_TARGET_PROC long double mc_legendre_q5l(const long double x)
{
//!# Legendre functions of the second kind, degree 5.
//!#
//!#    - Abramowitz, M. and Stegun, C. A. (Eds.). Legendre Functions. Ch. 8 in Handbook of Mathematical Functions with
//!#      Formulas, Graphs, and Mathematical Tables, 9th printing. New York: Dover, pp. 331-339, 1972.
//!#    - Arfken, G. Legendre Functions of the Second Kind, Qn(x). Mathematical Methods for Physicists, 3rd ed.
//!#      Orlando, FL: Academic Press, pp. 701-707, 1985.
//!#    - Binney, J. and Tremaine, S. Associated Legendre Functions. Appendix 5 in Galactic Dynamics. Princeton, NJ:
//!#      Princeton University Press, pp. 654-655, 1987.
//!#    - Morse, P. M. and Feshbach, H. Methods of Theoretical Physics, Part I. New York: McGraw-Hill, pp. 597-600, 1953.
//!#    - Snow, C. Hypergeometric and Legendre Functions with Applications to Integral Equations of Potential Theory.
//!#      Washington,  DC: U. S. Government Printing Office, 1952.
//!#    - Spanier, J. and Oldham, K. B. The Legendre Functions Pnu(x) and Qnu(x). Ch. 59 in An Atlas of Functions.
//!#      Washington, DC: Hemisphere, pp. 581-597, 1987.
//!#
	const long double p5 = mc_legendre_p5l(x);
	const long double q0 = mc_legendre_q0l(x);
	return p5 * q0 - (MCK_KL(MCK_63_8) * mc_raise4l(x)) + (MCK_KL(MCK_49_8) * mc_raise2l(x)) - MCK_KL(MCK_8_15);
}

#endif /* !MC_LEGENDRE_QN_H */

/* EOF */
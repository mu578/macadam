//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_math.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#ifndef MC_MATH_H
#define MC_MATH_H

#	include <macadam/details/math/mc_absmag.h>
#	include <macadam/details/math/mc_absrsqrt.h>
#	include <macadam/details/math/mc_abssqrt.h>
#	include <macadam/details/math/mc_acos.h>
#	include <macadam/details/math/mc_acosh.h>
#	include <macadam/details/math/mc_acospi.h>
#	include <macadam/details/math/mc_acot.h>
#	include <macadam/details/math/mc_acoth.h>
#	include <macadam/details/math/mc_acsc.h>
#	include <macadam/details/math/mc_acsch.h>
#	include <macadam/details/math/mc_asec.h>
#	include <macadam/details/math/mc_asech.h>
#	include <macadam/details/math/mc_asin.h>
#	include <macadam/details/math/mc_asinh.h>
#	include <macadam/details/math/mc_asinpi.h>
#	include <macadam/details/math/mc_atan.h>
#	include <macadam/details/math/mc_atan2.h>
#	include <macadam/details/math/mc_atan2pi.h>
#	include <macadam/details/math/mc_atanh.h>
#	include <macadam/details/math/mc_atanpi.h>
#	include <macadam/details/math/mc_bernoulli_b2n.h>
#	include <macadam/details/math/mc_bessel.h>
#	include <macadam/details/math/mc_besseli.h>
#	include <macadam/details/math/mc_besselj.h>
#	include <macadam/details/math/mc_bessely.h>
#	include <macadam/details/math/mc_beta.h>
#	include <macadam/details/math/mc_binomial.h>
#	include <macadam/details/math/mc_bisquare_psi.h>
#	include <macadam/details/math/mc_bisquare_rho.h>
#	include <macadam/details/math/mc_bisquare.h>
#	include <macadam/details/math/mc_biweight.h>
#	include <macadam/details/math/mc_cabs.h>
#	include <macadam/details/math/mc_cabs2.h>
#	include <macadam/details/math/mc_cacos.h>
#	include <macadam/details/math/mc_cacosh.h>
#	include <macadam/details/math/mc_cacot.h>
#	include <macadam/details/math/mc_cacoth.h>
#	include <macadam/details/math/mc_cacsc.h>
#	include <macadam/details/math/mc_cacsch.h>
#	include <macadam/details/math/mc_cadd.h>
#	include <macadam/details/math/mc_carg.h>
#	include <macadam/details/math/mc_casin.h>
#	include <macadam/details/math/mc_casinh.h>
#	include <macadam/details/math/mc_catan.h>
#	include <macadam/details/math/mc_catanh.h>
#	include <macadam/details/math/mc_cbrt.h>
#	include <macadam/details/math/mc_cbrti.h>
#	include <macadam/details/math/mc_cbrtu.h>
#	include <macadam/details/math/mc_ccos.h>
#	include <macadam/details/math/mc_ccosh.h>
#	include <macadam/details/math/mc_ccot.h>
#	include <macadam/details/math/mc_ccoth.h>
#	include <macadam/details/math/mc_ccsc.h>
#	include <macadam/details/math/mc_ccsch.h>
#	include <macadam/details/math/mc_cdiv.h>
#	include <macadam/details/math/mc_ceil.h>
#	include <macadam/details/math/mc_cexp.h>
#	include <macadam/details/math/mc_cexp2.h>
#	include <macadam/details/math/mc_cexpm1.h>
#	include <macadam/details/math/mc_chbevl.h>
#	include <macadam/details/math/mc_choose.h>
#	include <macadam/details/math/mc_cimag.h>
#	include <macadam/details/math/mc_ciseq.h>
#	include <macadam/details/math/mc_ciszero.h>
#	include <macadam/details/math/mc_clgamma.h>
#	include <macadam/details/math/mc_clog.h>
#	include <macadam/details/math/mc_clog10.h>
#	include <macadam/details/math/mc_clog1p.h>
#	include <macadam/details/math/mc_clog2.h>
#	include <macadam/details/math/mc_cmul.h>
#	include <macadam/details/math/mc_conj.h>
#	include <macadam/details/math/mc_copysign.h>
#	include <macadam/details/math/mc_cos.h>
#	include <macadam/details/math/mc_cosd.h>
#	include <macadam/details/math/mc_cosh.h>
#	include <macadam/details/math/mc_cospi.h>
#	include <macadam/details/math/mc_cot.h>
#	include <macadam/details/math/mc_coth.h>
#	include <macadam/details/math/mc_cpow.h>
#	include <macadam/details/math/mc_cpow10.h>
#	include <macadam/details/math/mc_cproj.h>
#	include <macadam/details/math/mc_creal.h>
#	include <macadam/details/math/mc_csc.h>
#	include <macadam/details/math/mc_csch.h>
#	include <macadam/details/math/mc_csin.h>
#	include <macadam/details/math/mc_csinh.h>
#	include <macadam/details/math/mc_csqrt.h>
#	include <macadam/details/math/mc_csub.h>
#	include <macadam/details/math/mc_ctan.h>
#	include <macadam/details/math/mc_ctanh.h>
#	include <macadam/details/math/mc_ctgamma.h>
#	include <macadam/details/math/mc_digamma.h>
#	include <macadam/details/math/mc_erf.h>
#	include <macadam/details/math/mc_erfc.h>
#	include <macadam/details/math/mc_eta.h>
#	include <macadam/details/math/mc_evalpoly.h>
#	include <macadam/details/math/mc_exp.h>
#	include <macadam/details/math/mc_exp10.h>
#	include <macadam/details/math/mc_exp10i.h>
#	include <macadam/details/math/mc_exp10m1.h>
#	include <macadam/details/math/mc_expi.h>
#	include <macadam/details/math/mc_expit.h>
#	include <macadam/details/math/mc_exp2.h>
#	include <macadam/details/math/mc_exp2i.h>
#	include <macadam/details/math/mc_exp2m1.h>
#	include <macadam/details/math/mc_expm1.h>
#	include <macadam/details/math/mc_fabs.h>
#	include <macadam/details/math/mc_factorial.h>
#	include <macadam/details/math/mc_fasttwosum.h>
#	include <macadam/details/math/mc_fdim.h>
#	include <macadam/details/math/mc_ffrac.h>
#	include <macadam/details/math/mc_fhrt.h>
#	include <macadam/details/math/mc_fiseq.h>
#	include <macadam/details/math/mc_fisint.h>
#	include <macadam/details/math/mc_fisnear.h>
#	include <macadam/details/math/mc_fisodd.h>
#	include <macadam/details/math/mc_fisval.h>
#	include <macadam/details/math/mc_fix.h>
#	include <macadam/details/math/mc_floor.h>
#	include <macadam/details/math/mc_fma.h>
#	include <macadam/details/math/mc_fmax.h>
#	include <macadam/details/math/mc_fmin.h>
#	include <macadam/details/math/mc_fmod.h>
#	include <macadam/details/math/mc_fmod2pi.h>
#	include <macadam/details/math/mc_fmodpi.h>
#	include <macadam/details/math/mc_frexp.h>
#	include <macadam/details/math/mc_gamma_p.h>
#	include <macadam/details/math/mc_gamma_q.h>
#	include <macadam/details/math/mc_gamma_r.h>
#	include <macadam/details/math/mc_gamma.h>
#	include <macadam/details/math/mc_gammaln_r.h>
#	include <macadam/details/math/mc_gammaln.h>
#	include <macadam/details/math/mc_gammaquo.h>
#	include <macadam/details/math/mc_gammasign.h>
#	include <macadam/details/math/mc_gcd.h>
#	include <macadam/details/math/mc_hermite_hen.h>
#	include <macadam/details/math/mc_hermite_hn.h>
#	include <macadam/details/math/mc_huber_loss.h>
#	include <macadam/details/math/mc_huber_psi.h>
#	include <macadam/details/math/mc_huber_rho.h>
#	include <macadam/details/math/mc_huber_weight.h>
#	include <macadam/details/math/mc_hurwitz_zeta.h>
#	include <macadam/details/math/mc_hypot.h>
#	include <macadam/details/math/mc_hypot3.h>
#	include <macadam/details/math/mc_hypot2.h>
#	include <macadam/details/math/mc_ibeta.h>
#	include <macadam/details/math/mc_iexp2.h>
#	include <macadam/details/math/mc_igamma.h>
#	include <macadam/details/math/mc_ilog2.h>
#	include <macadam/details/math/mc_ilogb.h>
#	include <macadam/details/math/mc_intge.h>
#	include <macadam/details/math/mc_inverf.h>
#	include <macadam/details/math/mc_invlogit.h>
#	include <macadam/details/math/mc_invprobit.h>
#	include <macadam/details/math/mc_isfinite.h>
#	include <macadam/details/math/mc_isinf.h>
#	include <macadam/details/math/mc_isnan.h>
#	include <macadam/details/math/mc_isnormal.h>
#	include <macadam/details/math/mc_itrunc.h>
#	include <macadam/details/math/mc_itrunc16.h>
#	include <macadam/details/math/mc_itrunc32.h>
#	include <macadam/details/math/mc_itrunc64.h>
#	include <macadam/details/math/mc_lbeta.h>
#	include <macadam/details/math/mc_lchoose.h>
#	include <macadam/details/math/mc_lcm.h>
#	include <macadam/details/math/mc_ldexp.h>
#	include <macadam/details/math/mc_legendre_pn.h>
#	include <macadam/details/math/mc_legendre_pnm.h>
#	include <macadam/details/math/mc_legendre_qn.h>
#	include <macadam/details/math/mc_lerp.h>
#	include <macadam/details/math/mc_lgamma.h>
#	include <macadam/details/math/mc_lgamma_r.h>
#	include <macadam/details/math/mc_llrint.h>
#	include <macadam/details/math/mc_llround.h>
#	include <macadam/details/math/mc_lmgamma.h>
#	include <macadam/details/math/mc_log.h>
#	include <macadam/details/math/mc_log10.h>
#	include <macadam/details/math/mc_log10p1.h>
#	include <macadam/details/math/mc_log1m.h>
#	include <macadam/details/math/mc_log1me.h>
#	include <macadam/details/math/mc_log1p.h>
#	include <macadam/details/math/mc_log1pe.h>
#	include <macadam/details/math/mc_log1pmx.h>
#	include <macadam/details/math/mc_log2.h>
#	include <macadam/details/math/mc_log2p1.h>
#	include <macadam/details/math/mc_logaddexp.h>
#	include <macadam/details/math/mc_logb.h>
#	include <macadam/details/math/mc_logbase.h>
#	include <macadam/details/math/mc_logcf.h>
#	include <macadam/details/math/mc_logdiffexp.h>
#	include <macadam/details/math/mc_logistic.h>
#	include <macadam/details/math/mc_logit.h>
#	include <macadam/details/math/mc_lognn.h>
#	include <macadam/details/math/mc_logodds.h>
#	include <macadam/details/math/mc_logp1.h>
#	include <macadam/details/math/mc_logradix.h>
#	include <macadam/details/math/mc_logsubexp.h>
#	include <macadam/details/math/mc_logx2pi.h>
#	include <macadam/details/math/mc_lrint.h>
#	include <macadam/details/math/mc_lround.h>
#	include <macadam/details/math/mc_maxmag.h>
#	include <macadam/details/math/mc_minmag.h>
#	include <macadam/details/math/mc_modf.h>
#	include <macadam/details/math/mc_nchoosek.h>
#	include <macadam/details/math/mc_nearbyint.h>
#	include <macadam/details/math/mc_nextafter.h>
#	include <macadam/details/math/mc_nexttoward.h>
#	include <macadam/details/math/mc_nthroot.h>
#	include <macadam/details/math/mc_polyroot2.h>
#	include <macadam/details/math/mc_polyroot3.h>
#	include <macadam/details/math/mc_pow.h>
#	include <macadam/details/math/mc_pow2.h>
#	include <macadam/details/math/mc_pow2i.h>
#	include <macadam/details/math/mc_pow10.h>
#	include <macadam/details/math/mc_powi.h>
#	include <macadam/details/math/mc_powm1.h>
#	include <macadam/details/math/mc_pown.h>
#	include <macadam/details/math/mc_probit.h>
#	include <macadam/details/math/mc_raise2.h>
#	include <macadam/details/math/mc_raise3.h>
#	include <macadam/details/math/mc_raise4.h>
#	include <macadam/details/math/mc_raise5.h>
#	include <macadam/details/math/mc_raise6.h>
#	include <macadam/details/math/mc_rem2pi_cw.h>
#	include <macadam/details/math/mc_rem90d.h>
#	include <macadam/details/math/mc_remainder.h>
#	include <macadam/details/math/mc_remint2.h>
#	include <macadam/details/math/mc_rempio2_cw.h>
#	include <macadam/details/math/mc_remquo.h>
#	include <macadam/details/math/mc_rgamma.h>
#	include <macadam/details/math/mc_riemann_zeta_n.h>
#	include <macadam/details/math/mc_riemann_zeta_p.h>
#	include <macadam/details/math/mc_riemann_zeta.h>
#	include <macadam/details/math/mc_rint.h>
#	include <macadam/details/math/mc_rootn.h>
#	include <macadam/details/math/mc_round.h>
#	include <macadam/details/math/mc_rsqr.h>
#	include <macadam/details/math/mc_rsqrt.h>
#	include <macadam/details/math/mc_scalb.h>
#	include <macadam/details/math/mc_scalbln.h>
#	include <macadam/details/math/mc_scalbn.h>
#	include <macadam/details/math/mc_sec.h>
#	include <macadam/details/math/mc_sech.h>
#	include <macadam/details/math/mc_sigmoid.h>
#	include <macadam/details/math/mc_signbit.h>
#	include <macadam/details/math/mc_sin.h>
#	include <macadam/details/math/mc_sinc.h>
#	include <macadam/details/math/mc_sincos.h>
#	include <macadam/details/math/mc_sincospi.h>
#	include <macadam/details/math/mc_sind.h>
#	include <macadam/details/math/mc_sinh.h>
#	include <macadam/details/math/mc_sinhcosh.h>
#	include <macadam/details/math/mc_sinpi.h>
#	include <macadam/details/math/mc_sqr.h>
#	include <macadam/details/math/mc_sqrt.h>
#	include <macadam/details/math/mc_sqrt1pm1.h>
#	include <macadam/details/math/mc_tan.h>
#	include <macadam/details/math/mc_tanh.h>
#	include <macadam/details/math/mc_tanpi.h>
#	include <macadam/details/math/mc_tgamma.h>
#	include <macadam/details/math/mc_trigamma.h>
#	include <macadam/details/math/mc_trunc.h>
#	include <macadam/details/math/mc_twoproduct.h>
#	include <macadam/details/math/mc_twosum.h>
#	include <macadam/details/math/mc_xchebevaln.h>
#	include <macadam/details/math/mc_xlog1px.h>
#	include <macadam/details/math/mc_xlog1py.h>
#	include <macadam/details/math/mc_xlogp1x.h>
#	include <macadam/details/math/mc_xlogp1y.h>
#	include <macadam/details/math/mc_xlogx.h>
#	include <macadam/details/math/mc_xlogy.h>
#	include <macadam/details/math/mc_xpolyevaln.h>
#	include <macadam/details/math/mc_zabs.h>
#	include <macadam/details/math/mc_zabs2.h>
#	include <macadam/details/math/mc_zadd.h>
#	include <macadam/details/math/mc_zacos.h>
#	include <macadam/details/math/mc_zacosh.h>
#	include <macadam/details/math/mc_zacot.h>
#	include <macadam/details/math/mc_zacoth.h>
#	include <macadam/details/math/mc_zacsc.h>
#	include <macadam/details/math/mc_zacsch.h>
#	include <macadam/details/math/mc_zatan.h>
#	include <macadam/details/math/mc_zatanh.h>
#	include <macadam/details/math/mc_zarg.h>
#	include <macadam/details/math/mc_zasin.h>
#	include <macadam/details/math/mc_zasinh.h>
#	include <macadam/details/math/mc_zconj.h>
#	include <macadam/details/math/mc_zcos.h>
#	include <macadam/details/math/mc_zcosh.h>
#	include <macadam/details/math/mc_zcot.h>
#	include <macadam/details/math/mc_zcoth.h>
#	include <macadam/details/math/mc_zcsc.h>
#	include <macadam/details/math/mc_zcsch.h>
#	include <macadam/details/math/mc_zdiv.h>
#	include <macadam/details/math/mc_zexp.h>
#	include <macadam/details/math/mc_zexp2.h>
#	include <macadam/details/math/mc_zexp10.h>
#	include <macadam/details/math/mc_zexpm1.h>
#	include <macadam/details/math/mc_zfadd.h>
#	include <macadam/details/math/mc_zfdiv.h>
#	include <macadam/details/math/mc_zfmul.h>
#	include <macadam/details/math/mc_zfpow.h>
#	include <macadam/details/math/mc_zfsub.h>
#	include <macadam/details/math/mc_zgamma.h>
#	include <macadam/details/math/mc_zgammaln.h>
#	include <macadam/details/math/mc_zinv.h>
#	include <macadam/details/math/mc_ziseq.h>
#	include <macadam/details/math/mc_ziszero.h>
#	include <macadam/details/math/mc_zlog.h>
#	include <macadam/details/math/mc_zlog1p.h>
#	include <macadam/details/math/mc_zlog2.h>
#	include <macadam/details/math/mc_zlog10.h>
#	include <macadam/details/math/mc_zmod.h>
#	include <macadam/details/math/mc_zmul.h>
#	include <macadam/details/math/mc_zneg.h>
#	include <macadam/details/math/mc_znorm.h>
#	include <macadam/details/math/mc_zpochhammer.h>
#	include <macadam/details/math/mc_zpolar.h>
#	include <macadam/details/math/mc_zpolyroot2.h>
#	include <macadam/details/math/mc_zpolyroot3.h>
#	include <macadam/details/math/mc_zpow.h>
#	include <macadam/details/math/mc_zpow10.h>
#	include <macadam/details/math/mc_zproj.h>
#	include <macadam/details/math/mc_zrecip.h>
#	include <macadam/details/math/mc_zsin.h>
#	include <macadam/details/math/mc_zsinh.h>
#	include <macadam/details/math/mc_zsqr.h>
#	include <macadam/details/math/mc_zsqrt.h>
#	include <macadam/details/math/mc_zsub.h>
#	include <macadam/details/math/mc_ztan.h>
#	include <macadam/details/math/mc_ztanh.h>

#endif /* !MC_MATH_H */

/* EOF */
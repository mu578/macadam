//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_lapack_lasr.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/lapack/mc_blas.h>
#include <macadam/details/math/mc_maxmag.h>

#ifndef MC_LAPACKE_LASR_H
#define MC_LAPACKE_LASR_H

#pragma mark - mc_lapack_slasr -

MC_TARGET_FUNC void mc_lapack_slasr(const char side, const char pivot, const char direct, const int m, const int n, const float * c, const float * s, float * a, const int lda)
{
	const float one = 1.0f, zero = 0.0f;

	int i, info, j;
	float ctemp, stemp, temp;

	info = 0;
	if (!(mc_blas_lsame(side, 'L') || mc_blas_lsame(side, 'R'))) {
		info = 1;
	} else if (!(mc_blas_lsame(pivot, 'V') || mc_blas_lsame(pivot, 'T') || mc_blas_lsame(pivot, 'B'))) {
		info = 2;
	} else if (!(mc_blas_lsame(direct, 'F') || mc_blas_lsame(direct, 'B'))) {
		info = 3;
	} else if (m < 0) {
		info = 4;
	} else if (n < 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, m)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("SLASR ", info);
		return;
	}

	if (m == 0 || n == 0) {
		return;
	}

	if (mc_blas_lsame(side, 'L')) {
		if (mc_blas_lsame(pivot, 'V')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 1; j <= (m - 1); ++j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                                   = mc_blas_matrix_at(a, lda, n, j + 1, i);
							mc_blas_matrix_at(a, lda, n, j + 1, i) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i)     = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, j, i);
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = (m - 1); j >= 1; --j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                                   = mc_blas_matrix_at(a, lda, n, j + 1, i);
							mc_blas_matrix_at(a, lda, n, j + 1, i) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i)     = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, j, i);
						}
					}
				}
			}
		} else if (mc_blas_lsame(pivot, 'T')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 2; j <= m; ++j) {
						ctemp = mc_blas_vector_at(c, j - 1);
						stemp = mc_blas_vector_at(s, j - 1);
						if (ctemp != one || stemp != zero) {
							for (i = 1; i <= n; ++i) {
								temp                               = mc_blas_matrix_at(a, lda, n, j, i);
								mc_blas_matrix_at(a, lda, n, j, i) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, 1, i);
								mc_blas_matrix_at(a, lda, n, 1, i) = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, 1, i);
							}
						}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = m; j >= 2; --j) {
					ctemp = mc_blas_vector_at(c, j - 1);
					stemp = mc_blas_vector_at(s, j - 1);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, 1, i);
							mc_blas_matrix_at(a, lda, n, 1, i) = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, 1, i);
						}
					}
				}
			}
		} else if (mc_blas_lsame(pivot, 'B')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 1; j <= (m - 1); ++j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i) = stemp * mc_blas_matrix_at(a, lda, n, m, i) + ctemp * temp;
							mc_blas_matrix_at(a, lda, n, m, i) = ctemp * mc_blas_matrix_at(a, lda, n, m, i) - stemp * temp;
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = m - 1; j >= 1; --j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i) = stemp * mc_blas_matrix_at(a, lda, n, m, i) + ctemp * temp;
							mc_blas_matrix_at(a, lda, n, m, i) = ctemp * mc_blas_matrix_at(a, lda, n, m, i) - stemp * temp;
						}
					}
				}
			}
		}
	} else if (mc_blas_lsame(side, 'R')) {
		if (mc_blas_lsame(pivot, 'V')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 1; j <= (n - 1); ++j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                                   = mc_blas_matrix_at(a, lda, n, i, j + 1);
							mc_blas_matrix_at(a, lda, n, i, j + 1) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j)     = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, i, j);
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = (n - 1); j >= 1; --j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                                   = mc_blas_matrix_at(a, lda, n, i, j + 1);
							mc_blas_matrix_at(a, lda, n, i, j + 1) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j)     = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, i, j);
						}
					}
				}
			}
		} else if (mc_blas_lsame(pivot, 'T')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 2; j <= n; ++j) {
					ctemp = mc_blas_vector_at(c, j - 1);
					stemp = mc_blas_vector_at(s, j - 1);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, i, 1);
							mc_blas_matrix_at(a, lda, n, i, 1) = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, i, 1);
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = n; j >= 2; --j) {
					ctemp = mc_blas_vector_at(c, j - 1);
					stemp = mc_blas_vector_at(s, j - 1);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, i, 1);
							mc_blas_matrix_at(a, lda, n, i, 1) = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, i, 1);
						}
					}
				}
			}
		} else if (mc_blas_lsame(pivot, 'B')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 1; j <= (n - 1); ++j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j) = stemp * mc_blas_matrix_at(a, lda, n, i, n) + ctemp * temp;
							mc_blas_matrix_at(a, lda, n, i, n) = ctemp * mc_blas_matrix_at(a, lda, n, i, n) - stemp * temp;
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = (n - 1); j >= 1; --j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j) = stemp * mc_blas_matrix_at(a, lda, n, i, n) + ctemp * temp;
							mc_blas_matrix_at(a, lda, n, i, n) = ctemp * mc_blas_matrix_at(a, lda, n, i, n) - stemp * temp;
						}
					}
				}
			}
		}
	}
}

#pragma mark - mc_lapack_dlasr -

MC_TARGET_FUNC void mc_lapack_dlasr(const char side, const char pivot, const char direct, const int m, const int n, const double * c, const double * s, double * a, const int lda)
{
	const double one = 1.0, zero = 0.0;

	int i, info, j;
	double ctemp, stemp, temp;

	info = 0;
	if (!(mc_blas_lsame(side, 'L') || mc_blas_lsame(side, 'R'))) {
		info = 1;
	} else if (!(mc_blas_lsame(pivot, 'V') || mc_blas_lsame(pivot, 'T') || mc_blas_lsame(pivot, 'B'))) {
		info = 2;
	} else if (!(mc_blas_lsame(direct, 'F') || mc_blas_lsame(direct, 'B'))) {
		info = 3;
	} else if (m < 0) {
		info = 4;
	} else if (n < 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, m)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("DLASR ", info);
		return;
	}

	if (m == 0 || n == 0) {
		return;
	}

	if (mc_blas_lsame(side, 'L')) {
		if (mc_blas_lsame(pivot, 'V')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 1; j <= (m - 1); ++j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                                   = mc_blas_matrix_at(a, lda, n, j + 1, i);
							mc_blas_matrix_at(a, lda, n, j + 1, i) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i)     = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, j, i);
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = (m - 1); j >= 1; --j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                                   = mc_blas_matrix_at(a, lda, n, j + 1, i);
							mc_blas_matrix_at(a, lda, n, j + 1, i) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i)     = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, j, i);
						}
					}
				}
			}
		} else if (mc_blas_lsame(pivot, 'T')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 2; j <= m; ++j) {
						ctemp = mc_blas_vector_at(c, j - 1);
						stemp = mc_blas_vector_at(s, j - 1);
						if (ctemp != one || stemp != zero) {
							for (i = 1; i <= n; ++i) {
								temp                               = mc_blas_matrix_at(a, lda, n, j, i);
								mc_blas_matrix_at(a, lda, n, j, i) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, 1, i);
								mc_blas_matrix_at(a, lda, n, 1, i) = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, 1, i);
							}
						}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = m; j >= 2; --j) {
					ctemp = mc_blas_vector_at(c, j - 1);
					stemp = mc_blas_vector_at(s, j - 1);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, 1, i);
							mc_blas_matrix_at(a, lda, n, 1, i) = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, 1, i);
						}
					}
				}
			}
		} else if (mc_blas_lsame(pivot, 'B')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 1; j <= (m - 1); ++j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i) = stemp * mc_blas_matrix_at(a, lda, n, m, i) + ctemp * temp;
							mc_blas_matrix_at(a, lda, n, m, i) = ctemp * mc_blas_matrix_at(a, lda, n, m, i) - stemp * temp;
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = m - 1; j >= 1; --j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i) = stemp * mc_blas_matrix_at(a, lda, n, m, i) + ctemp * temp;
							mc_blas_matrix_at(a, lda, n, m, i) = ctemp * mc_blas_matrix_at(a, lda, n, m, i) - stemp * temp;
						}
					}
				}
			}
		}
	} else if (mc_blas_lsame(side, 'R')) {
		if (mc_blas_lsame(pivot, 'V')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 1; j <= (n - 1); ++j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                                   = mc_blas_matrix_at(a, lda, n, i, j + 1);
							mc_blas_matrix_at(a, lda, n, i, j + 1) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j)     = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, i, j);
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = (n - 1); j >= 1; --j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                                   = mc_blas_matrix_at(a, lda, n, i, j + 1);
							mc_blas_matrix_at(a, lda, n, i, j + 1) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j)     = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, i, j);
						}
					}
				}
			}
		} else if (mc_blas_lsame(pivot, 'T')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 2; j <= n; ++j) {
					ctemp = mc_blas_vector_at(c, j - 1);
					stemp = mc_blas_vector_at(s, j - 1);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, i, 1);
							mc_blas_matrix_at(a, lda, n, i, 1) = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, i, 1);
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = n; j >= 2; --j) {
					ctemp = mc_blas_vector_at(c, j - 1);
					stemp = mc_blas_vector_at(s, j - 1);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, i, 1);
							mc_blas_matrix_at(a, lda, n, i, 1) = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, i, 1);
						}
					}
				}
			}
		} else if (mc_blas_lsame(pivot, 'B')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 1; j <= (n - 1); ++j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j) = stemp * mc_blas_matrix_at(a, lda, n, i, n) + ctemp * temp;
							mc_blas_matrix_at(a, lda, n, i, n) = ctemp * mc_blas_matrix_at(a, lda, n, i, n) - stemp * temp;
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = (n - 1); j >= 1; --j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j) = stemp * mc_blas_matrix_at(a, lda, n, i, n) + ctemp * temp;
							mc_blas_matrix_at(a, lda, n, i, n) = ctemp * mc_blas_matrix_at(a, lda, n, i, n) - stemp * temp;
						}
					}
				}
			}
		}
	}
}

#pragma mark - mc_lapack_llasr -

MC_TARGET_FUNC void mc_lapack_llasr(const char side, const char pivot, const char direct, const int m, const int n, const long double * c, const long double * s, long double * a, const int lda)
{
	const long double one = 1.0L, zero = 0.0L;

	int i, info, j;
	long double ctemp, stemp, temp;

	info = 0;
	if (!(mc_blas_lsame(side, 'L') || mc_blas_lsame(side, 'R'))) {
		info = 1;
	} else if (!(mc_blas_lsame(pivot, 'V') || mc_blas_lsame(pivot, 'T') || mc_blas_lsame(pivot, 'B'))) {
		info = 2;
	} else if (!(mc_blas_lsame(direct, 'F') || mc_blas_lsame(direct, 'B'))) {
		info = 3;
	} else if (m < 0) {
		info = 4;
	} else if (n < 0) {
		info = 5;
	} else if (lda < mc_maxmag(1, m)) {
		info = 9;
	}
	if (info != 0) {
		mc_blas_xerbla("LLASR ", info);
		return;
	}

	if (m == 0 || n == 0) {
		return;
	}

	if (mc_blas_lsame(side, 'L')) {
		if (mc_blas_lsame(pivot, 'V')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 1; j <= (m - 1); ++j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                                   = mc_blas_matrix_at(a, lda, n, j + 1, i);
							mc_blas_matrix_at(a, lda, n, j + 1, i) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i)     = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, j, i);
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = (m - 1); j >= 1; --j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                                   = mc_blas_matrix_at(a, lda, n, j + 1, i);
							mc_blas_matrix_at(a, lda, n, j + 1, i) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i)     = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, j, i);
						}
					}
				}
			}
		} else if (mc_blas_lsame(pivot, 'T')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 2; j <= m; ++j) {
						ctemp = mc_blas_vector_at(c, j - 1);
						stemp = mc_blas_vector_at(s, j - 1);
						if (ctemp != one || stemp != zero) {
							for (i = 1; i <= n; ++i) {
								temp                               = mc_blas_matrix_at(a, lda, n, j, i);
								mc_blas_matrix_at(a, lda, n, j, i) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, 1, i);
								mc_blas_matrix_at(a, lda, n, 1, i) = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, 1, i);
							}
						}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = m; j >= 2; --j) {
					ctemp = mc_blas_vector_at(c, j - 1);
					stemp = mc_blas_vector_at(s, j - 1);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, 1, i);
							mc_blas_matrix_at(a, lda, n, 1, i) = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, 1, i);
						}
					}
				}
			}
		} else if (mc_blas_lsame(pivot, 'B')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 1; j <= (m - 1); ++j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i) = stemp * mc_blas_matrix_at(a, lda, n, m, i) + ctemp * temp;
							mc_blas_matrix_at(a, lda, n, m, i) = ctemp * mc_blas_matrix_at(a, lda, n, m, i) - stemp * temp;
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = m - 1; j >= 1; --j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= n; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, j, i);
							mc_blas_matrix_at(a, lda, n, j, i) = stemp * mc_blas_matrix_at(a, lda, n, m, i) + ctemp * temp;
							mc_blas_matrix_at(a, lda, n, m, i) = ctemp * mc_blas_matrix_at(a, lda, n, m, i) - stemp * temp;
						}
					}
				}
			}
		}
	} else if (mc_blas_lsame(side, 'R')) {
		if (mc_blas_lsame(pivot, 'V')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 1; j <= (n - 1); ++j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                                   = mc_blas_matrix_at(a, lda, n, i, j + 1);
							mc_blas_matrix_at(a, lda, n, i, j + 1) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j)     = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, i, j);
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = (n - 1); j >= 1; --j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                                   = mc_blas_matrix_at(a, lda, n, i, j + 1);
							mc_blas_matrix_at(a, lda, n, i, j + 1) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j)     = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, i, j);
						}
					}
				}
			}
		} else if (mc_blas_lsame(pivot, 'T')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 2; j <= n; ++j) {
					ctemp = mc_blas_vector_at(c, j - 1);
					stemp = mc_blas_vector_at(s, j - 1);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, i, 1);
							mc_blas_matrix_at(a, lda, n, i, 1) = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, i, 1);
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = n; j >= 2; --j) {
					ctemp = mc_blas_vector_at(c, j - 1);
					stemp = mc_blas_vector_at(s, j - 1);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j) = ctemp * temp - stemp * mc_blas_matrix_at(a, lda, n, i, 1);
							mc_blas_matrix_at(a, lda, n, i, 1) = stemp * temp + ctemp * mc_blas_matrix_at(a, lda, n, i, 1);
						}
					}
				}
			}
		} else if (mc_blas_lsame(pivot, 'B')) {
			if (mc_blas_lsame(direct, 'F')) {
				for (j = 1; j <= (n - 1); ++j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j) = stemp * mc_blas_matrix_at(a, lda, n, i, n) + ctemp * temp;
							mc_blas_matrix_at(a, lda, n, i, n) = ctemp * mc_blas_matrix_at(a, lda, n, i, n) - stemp * temp;
						}
					}
				}
			} else if (mc_blas_lsame(direct, 'B')) {
				for (j = (n - 1); j >= 1; --j) {
					ctemp = mc_blas_vector_at(c, j);
					stemp = mc_blas_vector_at(s, j);
					if (ctemp != one || stemp != zero) {
						for (i = 1; i <= m; ++i) {
							temp                               = mc_blas_matrix_at(a, lda, n, i, j);
							mc_blas_matrix_at(a, lda, n, i, j) = stemp * mc_blas_matrix_at(a, lda, n, i, n) + ctemp * temp;
							mc_blas_matrix_at(a, lda, n, i, n) = ctemp * mc_blas_matrix_at(a, lda, n, i, n) - stemp * temp;
						}
					}
				}
			}
		}
	}
}

#endif /* !MC_LAPACKE_LASR_H */

/* EOF */
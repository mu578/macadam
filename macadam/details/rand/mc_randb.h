//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_randb.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/rand/mc_randi.h>

#ifndef MC_RANDB_H
#define MC_RANDB_H

#pragma mark - mc_randb_xtea -

MC_TARGET_PROC void mc_randb_xtea(const int r, const uint32_t key[4], uint32_t v[2])
{
//!# - r: number of rounds.
//!# - k: 128-bit key to encipher with.
//!# - v: 64-bit  block to encipher.
	int i      = 0;
	uint32_t s = 0U;
	uint32_t d = 0x9E3779B9U;
	for (; i < r; i++) {
		v[0] = v[0] + (((v[1] << 4) ^ (v[1] >> 5)) + v[1]) ^ (s + key[s & 3]);
		s    = s    + d;
		v[1] = v[1] + (((v[0] << 4) ^ (v[0] >> 5)) + v[0]) ^ (s + key[(s >> 11) & 3]);
	}
}

#pragma mark - mc_randb64_speck -

#	define mc_randb64_speck_ror(x, r) ((x >> r) | (x << (64U - r)))
#	define mc_randb64_speck_rol(x, r) ((x << r) | (x >> (64U - r)))
#	define mc_randb64_speck_r(x, y, k)       \
		(                                     \
			  x  = mc_randb64_speck_ror(x, 8U) \
			, x += y                           \
			, x ^= k                           \
			, y  = mc_randb64_speck_rol(y, 3U) \
			, y ^= x                           \
		)

MC_TARGET_PROC void mc_randb64_speck(const int r, const uint64_t key[2], uint64_t v[2])
{
//!# - r: number of rounds.
//!# - k: 128-bit key to encipher with.
//!# - v: 128-bit block to encipher.
	unsigned int i = 0;
	unsigned int n = r > 0 ? mc_cast_expr(unsigned int, r - 1) : 0;
	uint64_t y     = v[0];
	uint64_t x     = v[1];
	uint64_t b     = key[0];
	uint64_t a     = key[1];

	mc_randb64_speck_r(x, y, b);
	for (; i < n; i++) {
		mc_randb64_speck_r(a, b, i);
		mc_randb64_speck_r(x, y, b);
	}
	v[0] = y;
	v[1] = x;
}

#pragma mark - mc_randb_encipher -

MC_TARGET_PROC void mc_randb_encipher(const int n, const uint32_t key[4], uint32_t * v)
{
//!# By usage and in this very context, it renders described 64-bit
//!# cipher attacks unrelated to the matter being considered here.
//!# - n: v[n], `n` must be a multiple of 2 i.e `x` blocks of 2.
//!# - k: 128-bit key to encipher with.
//!# - v: blocks to encipher.
	int i        = 0;
	uint32_t * p = v;
	for (; i < n; i += 2) {
//!# Computing rounds, a minimum of 16 is performed.
		mc_randb_xtea(16 + 32 / n, key, p);
		p = p + 2;
	}
}

#pragma mark - mc_randb64_encipher -

MC_TARGET_PROC void mc_randb64_encipher(const int n, const uint64_t key[2], uint64_t * v)
{
//!# - n: v[n], `n` must be a multiple of 2 i.e `x` blocks of 2.
//!# - k: 128-bit key to encipher with.
//!# - v: blocks to encipher.
	int i        = 0;
	uint64_t * p = v;
	for (; i < n; i += 2) {
//!# Computing rounds, 32 is performed, set according to table.
		mc_randb64_speck(32, key, p);
		p = p + 2;
	}
}

#pragma mark - mc_randb_pack32be -

#	if MC_TARGET_C99 || MC_TARGET_CPP11
#	define mc_randb_pack32be(v, b)                    \
	mc_scope_begin                                    \
		(b)[0] = (((v) & UINT32_C(0xFF000000)) >> 24); \
		(b)[1] = (((v) & UINT32_C(0x00FF0000)) >> 16); \
		(b)[2] = (((v) & UINT32_C(0x0000FF00)) >>  8); \
		(b)[3] = (((v) & UINT32_C(0x000000FF)));       \
	mc_scope_end
#	else
#	define mc_randb_pack32be(v, b)                                        \
	mc_scope_begin                                                        \
		(b)[0] = (((v) & mc_cast_expr(const uint32_t, 0xFF000000)) >> 24); \
		(b)[1] = (((v) & mc_cast_expr(const uint32_t, 0x00FF0000)) >> 16); \
		(b)[2] = (((v) & mc_cast_expr(const uint32_t, 0x0000FF00)) >>  8); \
		(b)[3] = (((v) & mc_cast_expr(const uint32_t, 0x000000FF)));       \
	mc_scope_end
#	endif

#pragma mark - mc_randb64_pack64be -

#	if MC_TARGET_C99 || MC_TARGET_CPP11
#	define mc_randb64_pack64be(v, b)                          \
	mc_scope_begin                                            \
		(b)[0] = (((v) & UINT64_C(0xFF00000000000000)) >> 56); \
		(b)[1] = (((v) & UINT64_C(0x00FF000000000000)) >> 48); \
		(b)[2] = (((v) & UINT64_C(0x0000FF0000000000)) >> 40); \
		(b)[3] = (((v) & UINT64_C(0x000000FF00000000)) >> 32); \
		(b)[4] = (((v) & UINT64_C(0x00000000FF000000)) >> 24); \
		(b)[5] = (((v) & UINT64_C(0x0000000000FF0000)) >> 16); \
		(b)[6] = (((v) & UINT64_C(0x000000000000FF00)) >> 8);  \
		(b)[7] = (((v) & UINT64_C(0x00000000000000FF)));       \
	mc_scope_end
#	else
#	define mc_randb64_pack64be(v, b)                                              \
	mc_scope_begin                                                                \
		(b)[0] = (((v) & mc_cast_expr(const uint64_t, 0xFF00000000000000)) >> 56); \
		(b)[1] = (((v) & mc_cast_expr(const uint64_t, 0x00FF000000000000)) >> 48); \
		(b)[2] = (((v) & mc_cast_expr(const uint64_t, 0x0000FF0000000000)) >> 40); \
		(b)[3] = (((v) & mc_cast_expr(const uint64_t, 0x000000FF00000000)) >> 32); \
		(b)[4] = (((v) & mc_cast_expr(const uint64_t, 0x00000000FF000000)) >> 24); \
		(b)[5] = (((v) & mc_cast_expr(const uint64_t, 0x0000000000FF0000)) >> 16); \
		(b)[6] = (((v) & mc_cast_expr(const uint64_t, 0x000000000000FF00)) >> 8);  \
		(b)[7] = (((v) & mc_cast_expr(const uint64_t, 0x00000000000000FF)));       \
	mc_scope_end
#	endif

#pragma mark - mc_randb_randv -

MC_TARGET_PROC void mc_randb_randv(const int n, uint32_t * r)
{
	int i = 0;
	for (; i < n; i++) {
		r[i] = mc_randi32();
	}
}

#pragma mark - mc_randb64_randv -

MC_TARGET_PROC void mc_randb64_randv(const int n, uint64_t * r)
{
	int i = 0;
	for (; i < n; i++) {
		r[i] = mc_randi64();
	}
}

#pragma mark - mc_randb -

MC_TARGET_FUNC int mc_randb(const int n, void * dst)
{
//!# Generating `n` random bytes using a 64-bit block cipher and stores them into `dst`. It can be used
//!# as an entropy source for a CSPRNG (Cryptographically secure pseudorandom number generator).
	uint32_t block_m[12];
	uint8_t  block_b[32];

	uint8_t * p;
	uint32_t * block_k, * block_v;

	unsigned int offset;

	if (n > 0 && mc_nonnullptr(dst)) {
		p       = mc_cast(uint8_t *, dst);
		offset  = mc_cast(unsigned int, n);
		block_k = block_m + 0;
		block_v = block_m + 4;

		mc_randb_randv(4, block_k);
		while (offset > 32) {
			mc_randb_randv(8, block_v);
			mc_randb_encipher(8, block_k, block_v);
			mc_randb_pack32be(block_v[0], block_b + 0);
			mc_randb_pack32be(block_v[1], block_b + 4);
			mc_randb_pack32be(block_v[2], block_b + 8);
			mc_randb_pack32be(block_v[3], block_b + 12);
			mc_randb_pack32be(block_v[4], block_b + 16);
			mc_randb_pack32be(block_v[5], block_b + 20);
			mc_randb_pack32be(block_v[6], block_b + 24);
			mc_randb_pack32be(block_v[7], block_b + 28);
			mc_os_memcpy(p, block_b, 32);
			p      = p + 32;
			offset = offset - 32;
		}
		while (offset > 16) {
			mc_randb_randv(4, block_v);
			mc_randb_encipher(4, block_k, block_v);
			mc_randb_pack32be(block_v[0], block_b + 0);
			mc_randb_pack32be(block_v[1], block_b + 4);
			mc_randb_pack32be(block_v[2], block_b + 8);
			mc_randb_pack32be(block_v[3], block_b + 12);
			mc_os_memcpy(p, block_b, 16);
			p      = p + 16;
			offset = offset - 16;
		}
		while (offset > 8) {
			mc_randb_randv(2, block_v);
			mc_randb_encipher(2, block_k, block_v);
			mc_randb_pack32be(block_v[0], block_b + 0);
			mc_randb_pack32be(block_v[1], block_b + 4);
			mc_os_memcpy(p, block_b, 8);
			p      = p + 8;
			offset = offset - 8;
		}
		if (offset > 0) {
			mc_randb_randv(2, block_v);
			mc_randb_encipher(2U, block_k, block_v);
			mc_randb_pack32be(block_v[0], block_b + 0);
			mc_randb_pack32be(block_v[1], block_b + 4);
			mc_os_memcpy(p, block_b, offset);
		}
		return 0;
	}
	return -1;
}

#pragma mark - mc_randb64 -

MC_TARGET_FUNC int mc_randb64(const int n, void * dst)
{
//!# Generating `n` random bytes using a 128-bit block cipher and stores them into `dst`. It can be used
//!# as an entropy source for a CSPRNG (Cryptographically secure pseudorandom number generator).
	uint64_t block_m[10];
	uint8_t  block_b[64];

	uint8_t * p;
	uint64_t * block_k, * block_v;

	unsigned int offset;

	if (n > 0 && mc_nonnullptr(dst)) {
		p       = mc_cast(uint8_t *, dst);
		offset  = mc_cast(unsigned int, n);
		block_k = block_m + 0;
		block_v = block_m + 2;

		mc_randb64_randv(2, block_k);
		while (offset > 64) {
			mc_randb64_randv(8, block_v);
			mc_randb64_encipher(8, block_k, block_v);
			mc_randb64_pack64be(block_v[0], block_b + 0);
			mc_randb64_pack64be(block_v[1], block_b + 8);
			mc_randb64_pack64be(block_v[2], block_b + 16);
			mc_randb64_pack64be(block_v[3], block_b + 24);
			mc_randb64_pack64be(block_v[4], block_b + 32);
			mc_randb64_pack64be(block_v[5], block_b + 40);
			mc_randb64_pack64be(block_v[6], block_b + 48);
			mc_randb64_pack64be(block_v[7], block_b + 56);
			mc_os_memcpy(p, block_b, 64);
			p      = p + 64;
			offset = offset - 64;
		}
		while (offset > 32) {
			mc_randb64_randv(4, block_v);
			mc_randb64_encipher(4, block_k, block_v);
			mc_randb64_pack64be(block_v[0], block_b + 0);
			mc_randb64_pack64be(block_v[1], block_b + 8);
			mc_randb64_pack64be(block_v[2], block_b + 16);
			mc_randb64_pack64be(block_v[3], block_b + 24);
			mc_os_memcpy(p, block_b, 32);
			p      = p + 32;
			offset = offset - 32;
		}
		while (offset > 16) {
			mc_randb64_randv(2, block_v);
			mc_randb64_encipher(2, block_k, block_v);
			mc_randb64_pack64be(block_v[0], block_b + 0);
			mc_randb64_pack64be(block_v[1], block_b + 8);
			mc_os_memcpy(p, block_b, 16);
			p      = p + 16;
			offset = offset - 16;
		}
		if (offset > 0) {
			mc_randb64_randv(2, block_v);
			mc_randb64_encipher(2, block_k, block_v);
			mc_randb64_pack64be(block_v[0], block_b + 0);
			mc_randb64_pack64be(block_v[1], block_b + 8);
			mc_os_memcpy(p, block_b, offset);
		}
		return 0;
	}
	return -1;
}

#endif /* !MC_RANDB_H */

/* EOF */
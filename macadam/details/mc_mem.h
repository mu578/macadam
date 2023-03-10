//
// # -*- coding: utf-8, tab-width: 3 -*-

// mc_mem.h
//
// Copyright (C) 2019-2021 mu578. All rights reserved.
//

#include <macadam/details/mc_target.h>

#ifndef MC_MEM_H
#define MC_MEM_H

#pragma mark - mc_os_memset -

MC_TARGET_PROC void * mc_os_memset(void * b, int c, const  size_t len)
{
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	size_t max = sizeof(size_t) < sizeof(uint64_t) ?
		  mc_cast(size_t, UINT32_MAX)
		: mc_cast(size_t, UINT64_MAX)
	;
#	else
	size_t max = mc_cast(size_t, UINT32_MAX);
#	endif
	if (mc_nonnullptr(b) && len > 0 && len < max) {
#	if MC_TARGET_CPP98
		return ::memset(b, c, len);
#	else
		return memset(b, c, len);
#	endif
	}
	return MC_NULLPTR;
}

#pragma mark - mc_os_memzero -

MC_TARGET_PROC void * mc_os_memzero(void * b, const  size_t len)
{
	return mc_os_memset(b, 0, len);
}

#pragma mark - mc_os_memcpy -

MC_TARGET_PROC void * mc_os_memcpy(void * MC_TARGET_RESTRICT dest, const void * MC_TARGET_RESTRICT src, const  size_t len)
{
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	size_t max = sizeof(size_t) < sizeof(uint64_t) ?
		  mc_cast(size_t, UINT32_MAX)
		: mc_cast(size_t, UINT64_MAX)
	;
#	else
	size_t max = mc_cast(size_t, UINT32_MAX);
#	endif
	if (mc_nonnullptr(dest) && mc_nonnullptr(src) && (len > 0 && len < max)) {
#	if MC_TARGET_CPP98
		return ::memcpy(dest, src, len);
#	else
		return memcpy(dest, src, len);
#	endif
	}
	return MC_NULLPTR;
}

#pragma mark - mc_os_mempcpy -

MC_TARGET_PROC void * mc_os_mempcpy(void * MC_TARGET_RESTRICT dest, const void * MC_TARGET_RESTRICT src, const  size_t len)
{
#	if defined(__clang__)
#	if __has_builtin(__builtin_mempcpy)
	return __builtin_mempcpy(dest, src, len);
#	else
	return mc_cast_expr(uint8_t *, mc_os_mempcy(dest, src, len)) + len;
#	endif
#	elif defined(__GNUC__)
	return __builtin_mempcpy(dest, src, len);
#	else
	return mc_cast_expr(uint8_t *, mc_os_mempcy(dest, src, len)) + len;
#	endif
}

#pragma mark - mc_os_memmove -

MC_TARGET_PROC void * mc_os_memmove(void * MC_TARGET_RESTRICT dest, const void * MC_TARGET_RESTRICT src, const  size_t len)
{
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	size_t max = sizeof(size_t) < sizeof(uint64_t) ?
		  mc_cast(size_t, UINT32_MAX)
		: mc_cast(size_t, UINT64_MAX)
	;
#	else
	size_t max = mc_cast(size_t, UINT32_MAX);
#	endif
	if (mc_nonnullptr(dest) && mc_nonnullptr(src) && (len > 0 && len < max)) {
#	if MC_TARGET_CPP98
		return ::memmove(dest, src, len);
#	else
		return memmove(dest, src, len);
#	endif
	}
	return MC_NULLPTR;
}

#pragma mark - mc_os_memcmp -

MC_TARGET_PROC int mc_os_memcmp(const void * MC_TARGET_RESTRICT left, const void * MC_TARGET_RESTRICT right, const  size_t len)
{
#	if MC_TARGET_C99 || MC_TARGET_CPP11
	size_t max = sizeof(size_t) < sizeof(uint64_t) ?
		  mc_cast(size_t, UINT32_MAX)
		: mc_cast(size_t, UINT64_MAX)
	;
#	else
	size_t max = mc_cast(size_t, UINT32_MAX);
#	endif
	if (mc_nonnullptr(left) && mc_nonnullptr(right) && (len > 0 && len < max)) {
#	if MC_TARGET_CPP98
		return ::memcmp(left, right, len);
#	else
		return memcmp(left, right, len);
#	endif
	}
	return -1;
}

#pragma mark - mc_os_malloc -

MC_TARGET_PROC void * mc_os_malloc(size_t size)
{
#	if MC_TARGET_CPP98
	return ::malloc(size);
#	else
	return malloc(size);
#	endif
}

#pragma mark - mc_os_cmalloc -

MC_TARGET_PROC void * mc_os_cmalloc(size_t size)
{
	return mc_os_memzero(mc_os_malloc(size), size);
}

#pragma mark - mc_os_calloc -

MC_TARGET_PROC void * mc_os_calloc(size_t count, size_t size)
{
#	if MC_TARGET_CPP98
	return ::calloc(count, size);
#	else
	return calloc(count, size);
#	endif
}

#pragma mark - mc_os_realloc -

MC_TARGET_PROC void * mc_os_realloc(void * ptr, size_t size)
{
#	if MC_TARGET_CPP98
	return ::realloc(ptr, size);
#	else
	return realloc(ptr, size);
#	endif
}

#pragma mark - mc_os_free -

MC_TARGET_PROC void mc_os_free(void * ptr)
{
#	if MC_TARGET_CPP98
	::free(ptr);
#	else
	free(ptr);
#	endif
}

#pragma mark - mc_base_memset -

#	define mc_base_memset(src_addr, c, size) \
		mc_os_memset(src_addr, c, mc_cast_expr(size_t, size))

#pragma mark - mc_base_memzero -

#	define mc_base_memzero(src_addr, size) \
		mc_os_memzero(src_addr, mc_cast_expr(size_t, size))

#pragma mark - mc_base_memzero_type -

#	if MC_TARGET_CPP98 && !MC_TARGET_MEMSET_CC_OPTIM
#		define mc_base_memzero_type(type, n, x)                                                        \
		mc_scope_begin                                                                                 \
			if ((n) > 0) {                                                                              \
				::std::fill_n(x, n, mc_cast(type, 0));                                                   \
			}                                                                                           \
		mc_scope_end
#	elif MC_TARGET_CPP98 && MC_TARGET_MEMSET_CC_OPTIM
#		define mc_base_memzero_type(type, n, x)                                                        \
		mc_scope_begin                                                                                 \
			if ((n) > 0) {                                                                              \
				::memset(x, 0, mc_cast(size_t, n) * sizeof(type));                                       \
			}                                                                                           \
		mc_scope_end
#	elif MC_TARGET_C11 && !MC_TARGET_MEMSET_CC_OPTIM && defined(__STDC_LIB_EXT1__)
#		define mc_base_memzero_type(type, n, x)                                                        \
		mc_scope_begin                                                                                 \
			if ((n) > 0) {                                                                              \
				memset_s(x, mc_cast(size_t, n) * sizeof(type), 0, mc_cast(size_t, n) * sizeof(type));    \
			}                                                                                           \
		mc_scope_end
#	elif MC_TARGET_C99 && MC_TARGET_MEMSET_CC_OPTIM
#		define mc_base_memzero_type(type, n, x)                                                        \
		mc_scope_begin                                                                                 \
			if ((n) > 0) {                                                                              \
				memset(x, 0, mc_cast(size_t, n) * sizeof(type));                                         \
			}                                                                                           \
		mc_scope_end
#	else
#		define mc_base_memzero_type(type, n, x)                                                        \
		mc_scope_begin                                                                                 \
			size_t __mc_base_memzero_type_i = 0;                                                        \
			if ((n) > 0) {                                                                              \
				for (; __mc_base_memzero_type_i < mc_cast_expr(size_t, n); __mc_base_memzero_type_i++) { \
					x[__mc_base_memzero_type_i] = mc_cast_expr(type, 0);                                  \
				}                                                                                        \
			}                                                                                           \
		mc_scope_end
#	endif

#pragma mark - mc_base_memcpy -

#	define mc_base_memcpy(dest_addr, src_addr, size) \
		mc_os_memcpy(dest_addr, src_addr, mc_cast_expr(size_t, size))

#pragma mark - mc_base_memcpy_type -

#	if MC_TARGET_CPP98
#		define mc_base_memcpy_type(type, n, y, x)                                                             \
		mc_scope_begin                                                                                        \
			if ((n) > 0) {                                                                                     \
				::memcpy(y, x, mc_cast_expr(size_t, n) * sizeof(type));                                         \
			}                                                                                                  \
		mc_scope_end
#	elif MC_TARGET_C11 && defined(__STDC_LIB_EXT1__)
#		define mc_base_memcpy_type(type, n, y, x)                                                             \
		mc_scope_begin                                                                                        \
			if ((n) > 0) {                                                                                     \
				memcpy_s(y, mc_cast_expr(size_t, n) * sizeof(type), x, mc_cast_expr(size_t, n) * sizeof(type)); \
			}                                                                                                  \
		mc_scope_end
#	else
#		define mc_base_memcpy_type(type, n, y, x)                                                             \
		mc_scope_begin                                                                                        \
			size_t __mc_base_memcpy_type_i = 0;                                                                \
			if ((n) > 0) {                                                                                     \
				for (; __mc_base_memcpy_type_i < mc_cast_expr(size_t, n); __mc_base_memcpy_type_i++) {          \
					y[__mc_base_memcpy_type_i] = mc_cast_expr(type, x[__mc_base_memcpy_type_i]);                 \
				}                                                                                               \
			}                                                                                                  \
		mc_scope_end
#	endif

#pragma mark - mc_base_memmove -

#	define mc_base_memmove(dest_addr, src_addr, size) \
		mc_os_memmove(dest_addr, src_addr, mc_cast_expr(size_t, size))

#pragma mark - mc_base_memcmp -

#	define mc_base_memcmp(left_addr, right_addr, size) \
		mc_os_memcmp(left_addr, right_addr, mc_cast_expr(size_t, size))

#pragma mark - mc_base_sizeck -

#	define mc_base_sizeck(size) \
		((size) > 0 && mc_cast_expr(size_t, size) < mc_cast(size_t, MC_TARGET_ALLOCATOR_MAXSIZE) ? 0 : -1)

#pragma mark - mc_base_countck -

#	define mc_base_countck(item_type, count) \
		((count) > 0 && mc_cast_expr(size_t, count) * sizeof((item_type)) < mc_cast(size_t, MC_TARGET_ALLOCATOR_MAXSIZE) ? 0 : -1)

#pragma mark - mc_base_boundck -

#	define mc_base_boundck(offset, bound) \
	(offset) >= 0 && (bound) > 0 && AFG_CAST(size_t, (offset)) < AFG_CAST(size_t, (bound)))

#pragma mark - mc_base_alloc -

#	define mc_base_alloc(item_type, size) \
		mc_cast(item_type *, mc_os_malloc(mc_cast_expr(size_t, size)))

#	define mc_base_alloc_count(item_type, count) \
		 mc_cast(item_type *, mc_os_malloc(mc_cast_expr(size_t, count) * sizeof((item_type))))

#pragma mark - mc_base_calloc -

#	define mc_base_calloc(item_type, size) \
		mc_cast(item_type *, mc_os_cmalloc(mc_cast_expr(size_t, size)))

#	define mc_base_calloc_count(item_type, count) \
		mc_cast(item_type *, mc_os_calloc(mc_cast_expr(size_t, count), sizeof((item_type))))

#pragma mark - mc_base_realloc -

#	define mc_base_realloc(item_type, src_addr, newsize) \
		mc_cast(item_type *, mc_os_realloc(src_addr, mc_cast(size_t, newsize)))

#	define mc_base_realloc_count(item_type, src_addr, newcount) \
		mc_cast(item_type *, mc_os_realloc(src_addr, mc_cast_expr(size_t, newcount) * sizeof((item_type))))

#pragma mark - mc_alloc -

#	define mc_alloc(item_type, dest_addr, size)    \
	mc_scope_begin                                 \
		dest_addr = mc_base_alloc(item_type, size); \
	mc_scope_end

#	define mc_alloc_safe(item_type, dest_addr, size)                                        \
	mc_scope_begin                                                                          \
		dest_addr = mc_base_sizeck(size) == 0 ? mc_base_alloc(item_type, size) : MC_NULLPTR; \
	mc_scope_end

#	define mc_alloc_count(item_type, dest_addr, count)    \
	mc_scope_begin                                        \
		dest_addr = mc_base_alloc_count(item_type, count); \
	mc_scope_end

#	define mc_alloc_count_safe(item_type, dest_addr, count)                                                     \
	mc_scope_begin                                                                                              \
		dest_addr = mc_base_countck(item_type, count) == 0 ? mc_base_alloc_count(item_type, count) : MC_NULLPTR; \
	mc_scope_end

#pragma mark - mc_calloc -

#	define mc_calloc(item_type, dest_addr, size)    \
	mc_scope_begin                                  \
		dest_addr = mc_base_calloc(item_type, size); \
	mc_scope_end

#	define mc_calloc_safe(item_type, dest_addr, size)                                        \
	mc_scope_begin                                                                           \
		dest_addr = mc_base_sizeck(size) == 0 ? mc_base_calloc(item_type, size) : MC_NULLPTR; \
	mc_scope_end

#	define mc_calloc_count(item_type, dest_addr, count)    \
	mc_scope_begin                                         \
		dest_addr = mc_base_calloc_count(item_type, count); \
	mc_scope_end

#	define mc_calloc_count_safe(item_type, dest_addr, count)                                                     \
	mc_scope_begin                                                                                               \
		dest_addr = mc_base_countck(item_type, count) == 0 ? mc_base_calloc_count(item_type, count) : MC_NULLPTR; \
	mc_scope_end

#pragma mark - mc_realloc -

#	define mc_realloc(item_type, dest_addr, src_addr, newsize)    \
	mc_scope_begin                                                \
		dest_addr = mc_base_realloc(item_type, src_addr, newsize); \
	mc_scope_end

#	define mc_realloc_safe(item_type, dest_addr, src_addr, newsize)                        \
	mc_scope_begin                                                                         \
		if (0 == mc_base_sizeck(newsize)) {                                                 \
			if (MC_NULLPTR == (dest_addr = mc_base_realloc(item_type, src_addr, newsize))) { \
				if (mc_nonnullptr(src_addr)) {                                                \
					mc_os_free(src_addr);                                                      \
					src_addr = MC_NULLPTR;                                                     \
				}                                                                             \
			}                                                                                \
		} else {                                                                            \
			dest_addr = MC_NULLPTR;                                                          \
		}                                                                                   \
	mc_scope_end

#	define mc_realloc_count(item_type, dest_addr, src_addr, newcount)    \
	mc_scope_begin                                                       \
		dest_addr = mc_base_realloc_count(item_type, src_addr, newcount); \
	mc_scope_end

#	define mc_realloc_count_safe(item_type, dest_addr, src_addr, newcount)                        \
	mc_scope_begin                                                                                \
		if (0 == mc_base_countck(item_type, newcount)) {                                           \
			if (MC_NULLPTR == (dest_addr = mc_base_realloc_count(item_type, src_addr, newcount))) { \
				if (mc_nonnullptr(src_addr)) {                                                       \
					mc_os_free(src_addr);                                                             \
					src_addr = MC_NULLPTR;                                                            \
				}                                                                                    \
			}                                                                                       \
		} else {                                                                                   \
			dest_addr = MC_NULLPTR;                                                                 \
		}                                                                                          \
	mc_scope_end

#pragma mark - mc_dealloc -

#	define mc_dealloc(src_addr) \
		mc_os_free(src_addr)

#	define mc_dealloc_safe(src_addr)  \
	mc_scope_begin                    \
		if (mc_nonnullptr(src_addr)) { \
			mc_os_free(src_addr);       \
			src_addr = MC_NULLPTR;      \
		}                              \
	mc_scope_end

#endif /* !MC_MEM_H */

/* EOF */
/* SPDX-License-Identifier: MIT
 * Fast LZSA2 decoder. LZSA2 is Emmanuel Marty's format
 * (https://github.com/emmanuel-marty/lzsa); this decodes the framed stream
 * written by `lzsa -f 2` and is a drop-in for its decoder.
 *
 * Based on the LZSA2 decoder submitted to coreboot (change 95448, by Arthur
 * Heymans, MIT), with:
 *  - a fast path for the ~80% of commands without a literal or match-length
 *    escape;
 *  - 32-byte over-copies for matches (one AVX2 move, or four 8-byte ones);
 *  - every read and write bounds checked, so corrupt input returns SIZE_MAX.
 * Build it with -O2: -O3's inlining makes it 6-10% slower, -Os ~25%.
 */
#include <string.h>
#include "lzsa2.h"
#ifdef __AVX2__
#include <immintrin.h>
#endif

#ifndef __always_inline
#define __always_inline inline __attribute__((always_inline))
#endif
#define LZSA2_FRAME_MASK	0x7fffff
#define LZSA2_FRAME_UNCOMPRESSED 0x800000

struct nibble_state {
	uint8_t byte;
	int pending;
};

static __always_inline int lzsa2_nibble(const uint8_t **in, const uint8_t *in_end,
					struct nibble_state *ns)
{
	if (ns->pending == 0) {
		if (*in >= in_end)
			return -1;
		ns->byte = *(*in)++;
		ns->pending = 1;
		return ns->byte >> 4;
	}

	ns->pending = 0;
	return ns->byte & 0x0f;
}

static __always_inline int lzsa2_length(const uint8_t **in, const uint8_t *in_end,
					struct nibble_state *ns, uint32_t *len)
{
	int nib = lzsa2_nibble(in, in_end, ns);

	if (nib < 0)
		return -1;
	*len += nib;
	if (nib != 15)
		return 0;

	if (*in >= in_end)
		return -1;
	*len += *(*in)++;

	if (*len == 257) {
		if (*in + 2 > in_end)
			return -1;
		*len = (*in)[0] | (uint32_t)(*in)[1] << 8;
		*in += 2;
	} else if (*len == 256) {
		*len = 0;
	}

	return 0;
}

static __always_inline void lzsa2_copy8(void *dst, const void *src)
{
	uint64_t value;

	__builtin_memcpy(&value, src, sizeof(value));
	__builtin_memcpy(dst, &value, sizeof(value));
}

/* 32-byte copy, all loads before the stores: correct for offset >= length.
 * The caller checks that 32 bytes fit before the end of the output. */
static __always_inline void lzsa2_copy32(void *dst, const void *src)
{
#ifdef __AVX2__
	__m256i v = _mm256_loadu_si256((const __m256i *)src);
	_mm256_storeu_si256((__m256i *)dst, v);
#else
	uint64_t a, b, c, d;
	__builtin_memcpy(&a, (const uint8_t *)src + 0, 8);
	__builtin_memcpy(&b, (const uint8_t *)src + 8, 8);
	__builtin_memcpy(&c, (const uint8_t *)src + 16, 8);
	__builtin_memcpy(&d, (const uint8_t *)src + 24, 8);
	__builtin_memcpy((uint8_t *)dst + 0, &a, 8);
	__builtin_memcpy((uint8_t *)dst + 8, &b, 8);
	__builtin_memcpy((uint8_t *)dst + 16, &c, 8);
	__builtin_memcpy((uint8_t *)dst + 24, &d, 8);
#endif
}

static int lzsa2_block(const uint8_t **in, const uint8_t *in_end,
		       uint8_t **out, uint8_t *out_end, uint8_t *out_start)
{
	struct nibble_state ns = { 0 };
	uint32_t offset = 0;

	while (*in < in_end) {
		uint8_t token = *(*in)++;
		/* ~80% of commands have no literal or length escape: <= 2
		 * literals, a <= 8 byte match. With 64 bytes of room left their
		 * over-copies cannot leave the output, and the 2-byte literal read
		 * cannot leave the input because a 3-byte footer follows every
		 * block. */
		if (__builtin_expect((token & 0x18) != 0x18 && (token & 0x07) != 0x07 &&
				     out_end - *out >= 64, 1)) {
			uint32_t literals = (token >> 3) & 3;
			uint16_t lv;	/* 2-byte over-copy */
			__builtin_memcpy(&lv, *in, sizeof(lv));
			__builtin_memcpy(*out, &lv, sizeof(lv));
			*in += literals; *out += literals;
			if (*in >= in_end) break;
			switch (token & 0xc0) {
			case 0x00: {
				int v = lzsa2_nibble(in, in_end, &ns);
				if (v < 0) return -1;
				offset = ((((uint32_t)v << 1) | ((token >> 5) & 1)) ^ 0x1e) + 1;
				break;
			}
			case 0x40:
				if (*in >= in_end) return -1;
				offset = (((uint32_t)*(*in)++ | ((uint32_t)((token >> 5) & 1) << 8)) ^ 0x0ff) + 1;
				break;
			case 0x80: {
				int v = lzsa2_nibble(in, in_end, &ns);
				if (v < 0 || *in >= in_end) return -1;
				offset = ((uint32_t)*(*in)++ | ((uint32_t)v << 9) |
					  ((uint32_t)((token >> 5) & 1) << 8)) ^ 0x1eff;
				offset += 513;
				break;
			}
			default:
				if ((token & 0x20) == 0) {
					if (*in + 2 > in_end) return -1;
					offset = (((uint32_t)(*in)[0] << 8 | (*in)[1]) ^ 0xffff) + 1;
					*in += 2;
				}
				break;
			}
			uint32_t length = (token & 0x07) + 2;
			if (offset == 0 || offset > (uint32_t)(*out - out_start)) return -1;
			const uint8_t *src = *out - offset;
			uint8_t *dst = *out;
			if (offset >= length && out_end - dst >= 32) {
				lzsa2_copy32(dst, src);
			} else if (offset >= length && length <= 8 && out_end - dst >= 8) {
				lzsa2_copy8(dst, src);
			} else {
				size_t done = 0;
				if (offset >= 16) {
					while (done + 16 <= length) {
						lzsa2_copy8(dst + done, src + done);
						lzsa2_copy8(dst + done + 8, src + done + 8);
						done += 16;
					}
				}
				while (done < length) { dst[done] = src[done]; done++; }
			}
			*out += length;
			continue;
		}
		uint32_t literals = (token >> 3) & 3;

		if (literals != 3) {
			if (*in + 4 <= in_end && *out + 4 <= out_end) {
				uint32_t value;
				__builtin_memcpy(&value, *in, sizeof(value));
				__builtin_memcpy(*out, &value, sizeof(value));
			} else if (literals) {
				if (*in + literals > in_end || *out + literals > out_end)
					return -1;
				memcpy(*out, *in, literals);
			}
		} else {
			if (lzsa2_length(in, in_end, &ns, &literals) < 0)
				return -1;
			if (*in + literals > in_end || *out + literals > out_end)
				return -1;
			memcpy(*out, *in, literals);
		}
		*in += literals;
		*out += literals;

		if (*in >= in_end)
			break;

		switch (token & 0xc0) {
		case 0x00: {
			int v = lzsa2_nibble(in, in_end, &ns);
			if (v < 0)
				return -1;
			offset = ((((uint32_t)v << 1) | ((token >> 5) & 1)) ^ 0x1e) + 1;
			break;
		}
		case 0x40:
			if (*in >= in_end)
				return -1;
			offset = (((uint32_t)*(*in)++ | ((uint32_t)((token >> 5) & 1) << 8)) ^ 0x0ff) + 1;
			break;
		case 0x80: {
			int v = lzsa2_nibble(in, in_end, &ns);
			if (v < 0 || *in >= in_end)
				return -1;
			offset = ((uint32_t)*(*in)++ | ((uint32_t)v << 9) |
				  ((uint32_t)((token >> 5) & 1) << 8)) ^ 0x1eff;
			offset += 513;
			break;
		}
		default:
			if ((token & 0x20) == 0) {
				if (*in + 2 > in_end)
					return -1;
				offset = (((uint32_t)(*in)[0] << 8 | (*in)[1]) ^ 0xffff) + 1;
				*in += 2;
			} /* else: repeat previous offset */
			break;
		}

		uint32_t length = token & 0x07;
		if (length == 7) {
			length = 9;
			if (lzsa2_length(in, in_end, &ns, &length) < 0)
				return -1;
			if (length == 0)
				return -1;
		} else {
			length += 2;
		}

		if (offset == 0 || offset > (uint32_t)(*out - out_start))
			return -1;
		if (*out + length > out_end)
			return -1;

		const uint8_t *src = *out - offset;
		uint8_t *dst = *out;

		if (offset >= length && out_end - dst >= 32) {
			if (length <= 32) {
				lzsa2_copy32(dst, src);
			} else {
				size_t done = 0;
				while (done + 32 <= length && out_end - (dst + done) >= 32) {
					lzsa2_copy32(dst + done, src + done);
					done += 32;
				}
				if (done < length && out_end - (dst + done) >= 32) {
					lzsa2_copy32(dst + done, src + done);
				} else {
					while (done < length) { dst[done] = src[done]; done++; }
				}
			}
		} else if (offset >= length && length <= 8 && out_end - dst >= 8) {
			lzsa2_copy8(dst, src);
		} else {
			size_t done = 0;

			if (offset >= 16) {
				while (done + 16 <= length) {
					lzsa2_copy8(dst + done, src + done);
					lzsa2_copy8(dst + done + 8, src + done + 8);
					done += 16;
				}
			}
			while (done < length) {
				dst[done] = src[done];
				done++;
			}
		}
		*out += length;
	}

	return 0;
}

size_t lzsa2_decompress(const void *srcv, size_t srcn, void *dstv, size_t dstn)
{
	const uint8_t *in = srcv;
	const uint8_t *in_end = in + srcn;
	uint8_t *out = dstv;
	uint8_t *out_start = dstv;
	uint8_t *out_end = out + dstn;

	if (srcn < 4 || in[0] != 0x7b || in[1] != 0x9e || (in[2] & 0xe0) != 0x20)
		return SIZE_MAX;
	in += 3;

	while (in_end - in >= 3) {
		uint32_t frame = (uint32_t)in[0] | ((uint32_t)in[1] << 8) |
				 ((uint32_t)in[2] << 16);
		uint32_t size = frame & LZSA2_FRAME_MASK;

		in += 3;
		if (size == 0)	/* footer */
			break;
		if (size + 3 > (size_t)(in_end - in))	/* a footer must follow */
			return SIZE_MAX;

		if (frame & LZSA2_FRAME_UNCOMPRESSED) {
			if (size > (size_t)(out_end - out))
				return SIZE_MAX;
			memcpy(out, in, size);
			out += size;
			in += size;
		} else {
			if (lzsa2_block(&in, in + size, &out, out_end, out_start) < 0)
				return SIZE_MAX;
		}
	}

	return (size_t)(out - out_start);
}

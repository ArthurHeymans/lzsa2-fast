/* SPDX-License-Identifier: MIT */
#ifndef LZSA2_H
#define LZSA2_H

#include <stddef.h>
#include <stdint.h>

/*
 * Decompress an LZSA2 frame (as written by `lzsa -f 2`) from src[0..srcn)
 * into dst[0..dstn). Returns the decompressed size, or SIZE_MAX on a
 * malformed stream or when dst is too small. dst needs no slack: wide copies
 * are only used where they stay inside dst.
 */
size_t lzsa2_decompress(const void *src, size_t srcn, void *dst, size_t dstn);

#endif

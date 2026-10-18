/* SPDX-License-Identifier: MIT
 * lzsa2d in out [reps]: decompress an `lzsa -f 2` file; with reps, also time
 * the decoder (best of reps) and print bytes and MB/s. */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "lzsa2.h"

static uint8_t *slurp(const char *path, size_t *n)
{
	FILE *f = fopen(path, "rb");
	if (!f) { perror(path); exit(1); }
	fseek(f, 0, SEEK_END);
	*n = (size_t)ftell(f);
	fseek(f, 0, SEEK_SET);
	uint8_t *p = malloc(*n ? *n : 1);
	if (fread(p, 1, *n, f) != *n) { perror(path); exit(1); }
	fclose(f);
	return p;
}

int main(int argc, char **argv)
{
	if (argc < 3) { fprintf(stderr, "usage: lzsa2d in out [reps]\n"); return 1; }
	size_t n, cap = 64u << 20;
	uint8_t *in = slurp(argv[1], &n), *out = malloc(cap);
	size_t r = lzsa2_decompress(in, n, out, cap);
	if (r == SIZE_MAX) { fprintf(stderr, "%s: corrupt stream\n", argv[1]); return 2; }
	FILE *f = fopen(argv[2], "wb");
	if (!f || fwrite(out, 1, r, f) != r) { perror(argv[2]); return 1; }
	fclose(f);
	int reps = argc > 3 ? atoi(argv[3]) : 0;
	double best = 1e30;
	for (int i = 0; i < reps; i++) {
		struct timespec a, b;
		clock_gettime(CLOCK_MONOTONIC, &a);
		lzsa2_decompress(in, n, out, r);
		clock_gettime(CLOCK_MONOTONIC, &b);
		double t = (double)(b.tv_sec - a.tv_sec) + (double)(b.tv_nsec - a.tv_nsec) * 1e-9;
		if (t < best) best = t;
	}
	if (reps) printf("%zu -> %zu bytes, %.1f us, %.0f MB/s\n", n, r, best * 1e6, (double)r / best / 1e6);
	return 0;
}

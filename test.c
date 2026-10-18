/* SPDX-License-Identifier: MIT
 * test compressed original [...]: each compressed file must decode to its
 * original into an exactly sized buffer; then truncated and bit-flipped
 * copies must decode without touching memory outside the buffers (build
 * with -fsanitize=address,undefined to check that). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
	int fail = 0;
	srand(1);
	for (int a = 1; a + 1 < argc; a += 2) {
		size_t cn, on;
		uint8_t *c = slurp(argv[a], &cn), *o = slurp(argv[a + 1], &on);
		uint8_t *d = malloc(on ? on : 1);
		if (lzsa2_decompress(c, cn, d, on) != on || memcmp(d, o, on)) {
			printf("FAIL %s\n", argv[a]);
			fail = 1;
		}
		if (on && lzsa2_decompress(c, cn, d, on - 1) != SIZE_MAX) {
			printf("FAIL %s: short buffer accepted\n", argv[a]);
			fail = 1;
		}
		for (int i = 0; i < 2000; i++) {
			size_t n = i & 1 ? (size_t)rand() % cn : cn;
			uint8_t *x = malloc(n ? n : 1);
			memcpy(x, c, n);
			if (n > 3 && !(i & 1)) x[3 + (size_t)rand() % (n - 3)] ^= (uint8_t)(1u << (rand() & 7));
			lzsa2_decompress(x, n, d, on);
			free(x);
		}
		free(c); free(o); free(d);
	}
	puts(fail ? "FAILED" : "ok");
	return fail;
}

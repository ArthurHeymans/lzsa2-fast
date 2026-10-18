# lzsa2-fast

A fast, bounds-checked decoder for [LZSA2](https://github.com/emmanuel-marty/lzsa)
by Emmanuel Marty. It reads the framed stream that `lzsa -f 2` writes and can
replace upstream's decoder. On firmware it is about 10% faster, and corrupt
input cannot make it read or write outside its buffers.

It started as the LZSA2 decoder I submitted to coreboot (change 95448) and
was tuned while comparing firmware decompressors for coreboot.

## Numbers

Corpus: 28 coreboot stages (romstage, ramstage and postcar from 10 boards),
3,646,560 bytes, 1,795,384 bytes as LZSA2. Decode time is the sum of
per-stage medians over 7 passes on one core, built with gcc `-O2
-march=x86-64-v2`. Code size is `.text + .rodata` at x86-64 `-Os`.

| decoder | Intel N100 | Core i5-8350U | code |
|---|---:|---:|---:|
| upstream `lzsa_decompress_inmem` | 7.08 ms | 7.82 ms | 1.1 KB |
| lzsa2-fast | **6.39 ms** (−10%) | **6.86 ms** (−12%) | 1.9 KB |

coreboot's current LZSA2 decoder runs at upstream speed (6.52 vs 6.60 ms in a
second harness, where lzsa2-fast takes 6.07 ms).

## What makes it faster

- **Fast path:** about 80% of commands have no literal or match-length
  escape. When at least 64 bytes of output remain, those commands take a
  path whose over-copies need no further checks.
- **Wide copies:** matches use 32-byte over-copies, a single move with
  `-mavx2`.
- **Build with `-O2`:** for this loop `-O3` is 6-10% slower because of
  inlining, and `-Os` is about 15% slower.

All the other LZSA2 decoder ideas I measured were slower: a token lookup
table, branchless offset decoding and, on the N100, an if-chain in place of the
switch.

## Use

```c
#include "lzsa2.h"
size_t n = lzsa2_decompress(src, src_len, dst, dst_cap);  /* SIZE_MAX on error */
```

`dst` needs no slack: wide copies are used only where they stay inside it.
Compress with upstream: `lzsa -f 2 in out`.

```
make            # lzsa2d: lzsa2d in.lzsa2 out [reps] (reps: also time it)
make test       # test vectors + 2000 corruptions each under ASan/UBSan
```

## License

MIT, see [LICENSE](LICENSE). The LZSA2 format and its reference encoder are
Emmanuel Marty's ([zlib/CC0](https://github.com/emmanuel-marty/lzsa)).

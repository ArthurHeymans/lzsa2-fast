/* SPDX-License-Identifier: MIT
 * FHL: firmware Huffman LZ. A per-stream semi-static Huffman LZ77 meant to sit
 * between LZSA2 (tiny, zero RAM) and LZMA (dense, slow):
 *  - LSB-first bitstreams, single-level Huffman lookup tables (MAXL bits);
 *  - a joint literal-run/match-length command symbol (16x16);
 *  - LZMA-style distance slots and three rep distances;
 *  - literals coded with K tables selected by an encoder-clustered map of the
 *    literal two before (not the previous output byte, see below).
 *
 * Decoding runs in two passes. All literals are decoded first, into the tail
 * of the output buffer; they alternate between two bitstreams and each
 * stream's context is its own previous literal, so the two Huffman chains
 * run in parallel. The commands then copy literal runs from the tail and
 * matches from the output, with wide over-copies for the leading nfast
 * commands (counted by the encoder) for which those stay in bounds.
 *
 *   fhl c [-K k] [-p passes] [-d depth] in out   compress (and verify)
 *   fhl d in out                                 decompress
 *
 * Stream: u32 size | u8 K | u32 nlit | u32 nfast | u32 offsets of literal
 * streams A and B and of the distance stream | map | meta lengths | code
 * lengths | commands (code, literal-run and match-length extras) | the three
 * other streams. Every stream ends in 8 zero bytes of refill slack.
 * -DFHL_DECODER_ONLY leaves only fhl_decode(); -DFHL_LIB drops main().
 */
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fhl.h"

#ifndef MAXL
#define MAXL 10                     /* max Huffman code length = LUT bits; 10 keeps K tables in L1 */
#endif
#define NLIT 256
#define NCMD 256
#define NDIST 43
#define NMETA 15
#define MAXK 16                     /* literal tables in a stream */
#define MAXT 32                     /* encoder: K + Kf tables in the ANZ parse mode */
#define NLEN (MAXT * NLIT + NCMD + NDIST)
#define MIN_MATCH 2
#define HDR 25

static const uint32_t ll_base[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 14, 22, 38, 70, 134, 262};
static const uint8_t ll_bits[16] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 20};
static const uint32_t ml_base[16] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 24, 40, 72, 136, 264};
static const uint8_t ml_bits[16] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 16};

static unsigned log2k(unsigned k) { return k > 8 ? 4 : k > 4 ? 3 : k > 2 ? 2 : k > 1 ? 1 : 0; }

/* ------------------------------------------------------------------ */
/* decoder                                                             */
/* ------------------------------------------------------------------ */

struct br { uint64_t bb; unsigned bc; const uint8_t *p; };

static inline uint64_t rd64(const uint8_t *p) { uint64_t v; memcpy(&v, p, 8); return v; }
static inline uint32_t rd32le(const uint8_t *p) { return p[0] | p[1] << 8 | p[2] << 16 | (uint32_t)p[3] << 24; }
static inline void refill(struct br *b)
{
	b->bb |= rd64(b->p) << b->bc;
	b->p += (63 - b->bc) >> 3;
	b->bc |= 56;
}
static inline unsigned bits(struct br *b, unsigned n)
{
	unsigned v = (unsigned)(b->bb & ((1ull << n) - 1));
	b->bb >>= n; b->bc -= n;
	return v;
}
static inline unsigned sym(struct br *b, const uint16_t *lut, unsigned lb)
{
	unsigned e = lut[b->bb & ((1u << lb) - 1)];
	unsigned n = (e >> 8) & 15;
	b->bb >>= n; b->bc -= n;
	return e & 255;
}

/* Canonical codes, bit-reversed for an LSB-first reader, into rev[]. */
static void rcodes(const uint8_t *len, unsigned n, uint16_t *rev)
{
	unsigned cnt[MAXL + 2] = {0}, next[MAXL + 2], code = 0;
	for (unsigned s = 0; s < n; s++) cnt[len[s]]++;
	cnt[0] = 0;
	for (unsigned l = 1; l <= MAXL; l++) { code = (code + cnt[l - 1]) << 1; next[l] = code; }
	for (unsigned s = 0; s < n; s++) {
		unsigned l = len[s], x = l ? next[l]++ : 0;
		x = (x & 0x5555) << 1 | (x >> 1 & 0x5555);
		x = (x & 0x3333) << 2 | (x >> 2 & 0x3333);
		x = (x & 0x0f0f) << 4 | (x >> 4 & 0x0f0f);
		x = (x & 0x00ff) << 8 | (x >> 8 & 0x00ff);
		rev[s] = (uint16_t)(x >> (16 - l));   /* l == 0: unused */
	}
}
/* Fill a 2^lb LUT with ent(sym, len). A lone symbol has length 1 but only
 * one code: it gets both halves. */
#define FILL_LUT(lut, len, n, lb, ent) do {                                   \
	uint16_t rev_[256];                                                     \
	unsigned used_ = 0, one_ = 0;                                           \
	rcodes(len, n, rev_);                                                   \
	for (unsigned s_ = 0; s_ < (n); s_++) {                                \
		unsigned l_ = (len)[s_];                                        \
		if (!l_) continue;                                              \
		for (unsigned i_ = rev_[s_]; i_ < (1u << (lb)); i_ += 1u << l_) \
			(lut)[i_] = ent(s_, l_);                                \
		used_ += 1u << ((lb) - l_);                                     \
		one_ = s_;                                                      \
	}                                                                       \
	if (used_ == (1u << (lb)) >> 1)                                         \
		for (unsigned i_ = 0; i_ < (1u << (lb)); i_++) (lut)[i_] = ent(one_, 1); \
} while (0)

/* 16-bit entries sym | len << 8 | ctx[sym] << 12: literal tables carry the
 * context the decoded literal selects, so the next table needs no lookup. */
static void build_lut(uint16_t *lut, const uint8_t *len, unsigned n, unsigned lb, const uint8_t *ctx)
{
#define ENT16(s, l) (uint16_t)((s) | (l) << 8 | (ctx ? ctx[s] << 12 : 0))
	FILL_LUT(lut, len, n, lb, ENT16);
#undef ENT16
}

/* command: llx 5 | llb 9 | mlx 5 | mlb 9 (shifted by 4 for the length) */
static uint32_t cmd_payload(unsigned c)
{
	unsigned ll = c >> 4, ml = c & 15;
	return ll_bits[ll] | ll_base[ll] << 5 | (uint32_t)ml_bits[ml] << 14 | (uint32_t)ml_base[ml] << 19;
}
/* distance: rep 2 (3 = new) | nb 5 | base 21 */
static uint32_t dist_payload(unsigned ds)
{
	if (ds < 3) return ds;
	unsigned slot = ds - 3;
	if (slot < 4) return 3 | (uint32_t)(slot + 1) << 7;
	unsigned nb = (slot >> 1) - 1;
	return 3 | nb << 2 | (((2u | (slot & 1)) << nb) + 1) << 7;
}

/* 32-bit LUT for command/distance symbols: len in bits 0-3, payload above */
static void build_lut32(uint32_t *lut, const uint8_t *len, unsigned n, uint32_t (*payload)(unsigned sym))
{
	uint32_t pay[256];
	for (unsigned s = 0; s < n; s++) pay[s] = payload(s) << 4;
#define ENT32(s, l) ((l) | pay[s])
	FILL_LUT(lut, len, n, MAXL, ENT32);
#undef ENT32
}

struct dtab {
	uint16_t lit[MAXK][1 << MAXL];
	uint32_t cmd[1 << MAXL], dist[1 << MAXL];
	uint64_t msk[33];
	uint8_t len[MAXK * NLIT + NCMD + NDIST];
	uint32_t size;
	unsigned lit0;          /* table of the initial context */
};

_Static_assert(sizeof(struct dtab) <= FHL_WORK_SIZE, "FHL_WORK_SIZE too small");

/* Header and tables; returns the bit position of the command stream. */
__attribute__((noinline)) static uint64_t fhl_setup(const uint8_t *src, struct dtab *t)
{
	t->size = rd32le(src);
	unsigned K = src[4], kb = log2k(K);
	struct br b = { 0, 0, src + HDR };
	uint8_t map[256];
	for (unsigned i = 0; i < 256; i++) {
		refill(&b);
		map[i] = (uint8_t)bits(&b, kb);
	}
	uint8_t mlen[NMETA];
	uint16_t mlut[1 << 7];
	refill(&b);
	for (unsigned i = 0; i < NMETA; i++) mlen[i] = (uint8_t)bits(&b, 3);
	build_lut(mlut, mlen, NMETA, 7, NULL);
	unsigned total = K * NLIT + NCMD + NDIST, i = 0, prev = 0;
	uint8_t *len = t->len;
	while (i < total) {
		refill(&b);
		unsigned s = sym(&b, mlut, 7), r;
		if (s < 12) { len[i++] = (uint8_t)(prev = s); continue; }
		if (s == 12) { r = 3 + bits(&b, 2); while (r--) len[i++] = (uint8_t)prev; continue; }
		r = s == 13 ? 3 + bits(&b, 3) : 11 + bits(&b, 7);
		while (r--) len[i++] = 0;
		prev = 0;
	}
	for (unsigned k = 0; k < K; k++) build_lut(t->lit[k], len + k * NLIT, NLIT, MAXL, map);
	build_lut32(t->cmd, len + K * NLIT, NCMD, cmd_payload);
	build_lut32(t->dist, len + K * NLIT + NCMD, NDIST, dist_payload);
	for (unsigned k = 0; k <= 32; k++) t->msk[k] = k ? ~0ull >> (64 - k) : 0;
	t->lit0 = map[0];
	return (uint64_t)(b.p - src) * 8 - b.bc;
}

/* gcc's SLP pass turns the rep array below into xmm inserts: slower */
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("no-tree-slp-vectorize")))
#endif
__attribute__((noinline)) size_t fhl_decode(const void *srcv, size_t n, void *dstv, size_t cap, void *work)
{
	const uint8_t *src = srcv;
	uint8_t *dst = dstv;
	struct dtab *t = work;
	if (n < HDR || src[4] == 0 || src[4] > MAXK) return 0;
	uint64_t cpos = fhl_setup(src, t);
	uint32_t size = t->size, nlit = rd32le(src + 5), nfast = rd32le(src + 9);
	if (size > cap || nlit > size) return 0;
	uint8_t *o = dst, *const end = dst + size;
	const uint8_t *lr = end - nlit;

	/* pass 1: literal 2i from stream A, 2i+1 from stream B, into the tail */
	{
		const uint8_t *pa = src + rd32le(src + 13), *pb = src + rd32le(src + 17);
		const uint16_t *const lits = t->lit[0];
		const uint16_t *ta = t->lit[t->lit0], *tb = ta;
		uint64_t ba = 0, bb = 0;
		unsigned ca = 0, cb = 0;
		uint8_t *q = end - nlit;
		uint32_t k = nlit;
#define RF(b, c, p) do { b |= rd64(p) << c; p += (63 - c) >> 3; c |= 56; } while (0)
#define LIT(b, c, tt) do {                                                     \
		unsigned e_ = tt[b & ((1u << MAXL) - 1)], k_ = (e_ >> 8) & 15; \
		b >>= k_; c -= k_;                                            \
		*q++ = (uint8_t)e_;                                           \
		tt = lits + ((size_t)(e_ >> 12) << MAXL);                     \
	} while (0)
		while (k >= 10) {          /* 5 codes of <= MAXL bits per refill */
			RF(ba, ca, pa); RF(bb, cb, pb);
			LIT(ba, ca, ta); LIT(bb, cb, tb); LIT(ba, ca, ta); LIT(bb, cb, tb);
			LIT(ba, ca, ta); LIT(bb, cb, tb); LIT(ba, ca, ta); LIT(bb, cb, tb);
			LIT(ba, ca, ta); LIT(bb, cb, tb);
			k -= 10;
		}
		RF(ba, ca, pa); RF(bb, cb, pb);
		while (k >= 2) { LIT(ba, ca, ta); LIT(bb, cb, tb); k -= 2; }
		if (k) LIT(ba, ca, ta);
#undef RF
#undef LIT
	}

	/* pass 2: commands. A 64-bit peek covers a command code and both length
	 * extras; distances come from their own stream. */
	uint64_t dpos = (uint64_t)rd32le(src + 21) * 8;
	uint32_t r0 = 1, r1 = 1, r2 = 1;
	for (;;) {
		uint64_t pk = rd64(src + (cpos >> 3)) >> (cpos & 7);
		uint32_t ce = t->cmd[pk & ((1u << MAXL) - 1)];
		unsigned l0 = ce & 15, llx = (ce >> 4) & 31, mlx = (ce >> 18) & 31;
		uint32_t L = ((ce >> 9) & 511) + (uint32_t)((pk >> l0) & t->msk[llx]);
		uint32_t M = (ce >> 23) + (uint32_t)((pk >> (l0 + llx)) & t->msk[mlx]);
		cpos += l0 + llx + mlx;
		if (__builtin_expect(nfast != 0, 1)) {
			memcpy(o, lr, 16);
			if (L > 16) memmove(o + 16, lr + 16, L - 16);
			o += L; lr += L;
		} else {
			while (L--) *o++ = *lr++;
			if (o >= end) break;
		}
		uint64_t dk = rd64(src + (dpos >> 3)) >> (dpos & 7);
		uint32_t de = t->dist[dk & ((1u << MAXL) - 1)];
		unsigned d0 = de & 15, rep = (de >> 4) & 3, nb = (de >> 6) & 31;
		uint32_t far = (de >> 11) + (uint32_t)((dk >> d0) & t->msk[nb]);
		dpos += d0 + nb;
		uint32_t rr[4] = { r0, r1, r2, far };  /* indexed: a branch chain mispredicts */
		uint32_t d = rr[rep];
		r2 = rep >= 2 ? r1 : r2;
		r1 = rep >= 1 ? r0 : r1;
		r0 = d;
		const uint8_t *s = o - d;
		uint8_t *e = o + M;
		if (__builtin_expect(nfast != 0 && d >= 8, 1)) {
			/* two 8-byte halves: valid for any d >= 8, covers M <= 16 */
			memcpy(o, s, 8);
			memcpy(o + 8, s + 8, 8);
			if (__builtin_expect(M > 16, 0)) {
				o += 16; s += 16;
				do { memcpy(o, s, 8); o += 8; s += 8; } while (o < e);
			}
		} else {
			while (o < e) *o++ = *s++;
		}
		o = e;
		nfast -= nfast != 0;
	}
	return (size_t)(o - dst);
}

#ifndef FHL_DECODER_ONLY
/* ------------------------------------------------------------------ */
/* encoder                                                             */
/* ------------------------------------------------------------------ */

struct bw { uint8_t *p; uint64_t acc; unsigned cnt; };
static void put(struct bw *w, uint32_t v, unsigned n)
{
	w->acc |= (uint64_t)v << w->cnt;
	w->cnt += n;
	while (w->cnt >= 8) { *w->p++ = (uint8_t)w->acc; w->acc >>= 8; w->cnt -= 8; }
}

/* Huffman code lengths limited to maxl (frequency halving until it fits). */
static void huff_len(const uint32_t *freq, unsigned n, unsigned maxl, uint8_t *len)
{
	uint64_t f[512];
	int par[512], act[512];
	for (unsigned scale = 0;; scale++) {
		int nn = 0, na = 0;
		memset(len, 0, n);
		for (unsigned s = 0; s < n; s++)
			if (freq[s]) { f[nn] = scale ? (freq[s] >> scale) + 1 : freq[s]; par[nn] = -1; act[na++] = nn++; }
		if (nn == 0) return;
		if (nn == 1) { for (unsigned s = 0; s < n; s++) if (freq[s]) len[s] = 1; return; }
		while (na > 1) {
			int a = 0, b2 = 1;
			if (f[act[b2]] < f[act[a]]) { a = 1; b2 = 0; }
			for (int i = 2; i < na; i++) {
				if (f[act[i]] < f[act[a]]) { b2 = a; a = i; }
				else if (f[act[i]] < f[act[b2]]) b2 = i;
			}
			f[nn] = f[act[a]] + f[act[b2]]; par[nn] = -1;
			par[act[a]] = nn; par[act[b2]] = nn;
			int hi = a > b2 ? a : b2, lo = a > b2 ? b2 : a;
			act[hi] = act[--na]; act[lo] = nn++;
		}
		unsigned mx = 0, k = 0;
		for (unsigned s = 0; s < n; s++) {
			if (!freq[s]) continue;
			unsigned d = 0;
			for (int x = k; par[x] >= 0; x = par[x]) d++;
			len[s] = (uint8_t)d; if (d > mx) mx = d;
			k++;
		}
		if (mx <= maxl) return;
	}
}

static void codes(const uint8_t *len, unsigned n, unsigned maxl, uint32_t *code)
{
	unsigned cnt[16] = {0}, next[16], c = 0;
	for (unsigned s = 0; s < n; s++) cnt[len[s]]++;
	cnt[0] = 0;
	for (unsigned l = 1; l <= maxl; l++) { c = (c + cnt[l - 1]) << 1; next[l] = c; }
	for (unsigned s = 0; s < n; s++) {
		unsigned l = len[s], cc = l ? next[l]++ : 0, r = 0;
		for (unsigned i = 0; i < l; i++) r |= ((cc >> i) & 1) << (l - 1 - i);
		code[s] = r;
	}
}

static unsigned dist_slot(uint32_t v)
{
	if (v < 4) return v;
	unsigned n = 31 - __builtin_clz(v);
	return 2 * n + ((v >> (n - 1)) & 1);
}
static unsigned ll_code(uint32_t L) { int c = 15; while (L < ll_base[c]) c--; return c; }
static unsigned ml_code(uint32_t M) { int c = 15; while (c && M < ml_base[c]) c--; return c; }

/* ds: 0..2 rep, else 3 + distance slot; d the actual distance */
#define seq fhl_seq

struct model {
	unsigned K, Kf;
	uint8_t map[256], mapf[256];
};

struct cost {                         /* 1/16 bit */
	uint32_t lit[MAXT][256], cmd[256], dist[NDIST];
};

static void lens_to_cost(const uint8_t *len, unsigned n, uint32_t *c, unsigned miss)
{
	for (unsigned s = 0; s < n; s++) c[s] = len[s] ? len[s] * 16 : miss * 16;
}

struct cand { uint32_t len, dist; };
struct mf { struct cand *c; uint32_t *off; };

static uint32_t mlen(const uint8_t *a, const uint8_t *b, uint32_t max)
{
	uint32_t l = 0;
	while (l + 8 <= max) {
		uint64_t x = rd64(a + l) ^ rd64(b + l);
		if (x) return l + (__builtin_ctzll(x) >> 3);
		l += 8;
	}
	while (l < max && a[l] == b[l]) l++;
	return l;
}

#define MAXM (264 + 65535)
#define NICE 1024
#define MAXD (1u << 20)
static void find_matches(const uint8_t *in, uint32_t n, unsigned depth, struct mf *m)
{
	enum { HB = 18 };
	uint32_t *head = malloc(sizeof(uint32_t) << HB), *prv = malloc(4 * (size_t)n + 4);
	memset(head, 0xff, sizeof(uint32_t) << HB);
	size_t cap = (size_t)n * 4 + 1024, k = 0;
	m->c = malloc(cap * sizeof(struct cand));
	m->off = malloc(4 * ((size_t)n + 1));
	uint32_t run_end = 0, run_d = 0;   /* a match >= NICE covers i..run_end */
	for (uint32_t i = 0; i < n; i++) {
		m->off[i] = (uint32_t)k;
		if (i + 3 > n) continue;
		uint32_t h = ((in[i] | in[i + 1] << 8 | in[i + 2] << 16) * 2654435761u) >> (32 - HB);
		uint32_t j = head[h], best = 2, max = n - i < MAXM ? n - i : MAXM;
		head[h] = i; prv[i] = j;
		if (run_end >= i + NICE) {
			/* inside a long match: its continuation is the only candidate
			 * (the parser jumps over it anyway); no chain walk, no compare */
			if (k + 1 >= cap) { cap *= 2; m->c = realloc(m->c, cap * sizeof(struct cand)); }
			m->c[k].len = run_end - i < max ? run_end - i : max; m->c[k].dist = run_d; k++;
			continue;
		}
		for (unsigned dd = 0; dd < depth && j != 0xffffffff; dd++, j = prv[j]) {
			if (i - j > MAXD) break;
			if (in[j + best] != in[i + best]) continue;
			uint32_t l = mlen(in + i, in + j, max);
			if (l > best) {
				if (k + 1 >= cap) { cap *= 2; m->c = realloc(m->c, cap * sizeof(struct cand)); }
				m->c[k].len = l; m->c[k].dist = i - j; k++;
				best = l;
				if (l >= NICE) { run_end = i + l; run_d = i - j; }
				if (l >= NICE || l == max) break;
			}
		}
	}
	m->off[n] = (uint32_t)k;
	free(head); free(prv);
}

/* fhl_parse() for ANZ prices literals with ANZ's context instead: the
 * previous output byte, and with Kf > 0 the match byte for the first literal
 * after a match */
static int legacy;

/* literal table index at position pos; plit holds the last two literals.
 * FHL's context is the literal two back: each literal stream then only
 * depends on itself. */
static inline unsigned lctx(const struct model *md, const uint8_t *in, uint32_t pos, int first,
			    uint32_t rep0, unsigned plit)
{
	if (!legacy) return md->map[plit >> 8];
	if (first && md->Kf) return md->K + md->mapf[in[pos - rep0]];
	return md->map[pos ? in[pos - 1] : 0];
}

struct node { uint32_t cost, from, M, ds, d, lr, rep[3]; uint16_t plit; uint8_t first; };

static size_t parse(const uint8_t *in, uint32_t n, const struct mf *m, const struct model *md,
		    const struct cost *c, struct seq *out)
{
	struct node *nd = malloc(sizeof(*nd) * ((size_t)n + 1));
	for (uint32_t i = 0; i <= n; i++) nd[i].cost = 0xffffffff;
	nd[0] = (struct node){ .cost = 0, .rep = { 1, 1, 1 } };
	for (uint32_t i = 0; i < n; i++) {
		struct node *a = &nd[i];
		uint32_t base = a->cost;
		unsigned llc = ll_code(a->lr);
		uint32_t lc = c->lit[lctx(md, in, i, a->first && a->lr == 0, a->rep[0], a->plit)][in[i]];
		if (base + lc < nd[i + 1].cost) {
			struct node *z = &nd[i + 1];
			z->cost = base + lc; z->from = i; z->M = 0; z->lr = a->lr + 1; z->first = a->first;
			z->plit = (uint16_t)(a->plit << 8 | in[i]);
			memcpy(z->rep, a->rep, 12);
		}
		uint32_t max = n - i < MAXM ? n - i : MAXM, longest = 0;
		uint32_t cmdbase = base + ll_bits[llc] * 16;
		for (unsigned r = 0; r < 3; r++) {
			uint32_t d = a->rep[r];
			if (d > i) continue;
			uint32_t l = mlen(in + i, in + i - d, max);
			longest = l > longest ? l : longest;
			if (l < MIN_MATCH) continue;
			for (uint32_t x = MIN_MATCH; x <= l; x++) {
				if (x > 64 && x != l) continue;
				unsigned mc = ml_code(x);
				uint32_t t = cmdbase + c->cmd[llc * 16 + mc] + ml_bits[mc] * 16 + c->dist[r];
				struct node *z = &nd[i + x];
				if (t < z->cost) {
					z->cost = t; z->from = i; z->M = x; z->ds = r; z->d = d; z->lr = 0; z->first = 1;
					z->plit = a->plit;
					if (r == 0) memcpy(z->rep, a->rep, 12);
					else if (r == 1) { z->rep[0] = d; z->rep[1] = a->rep[0]; z->rep[2] = a->rep[2]; }
					else { z->rep[0] = d; z->rep[1] = a->rep[0]; z->rep[2] = a->rep[1]; }
				}
			}
		}
		uint32_t lo = MIN_MATCH + 1;
		for (uint32_t k = m->off[i]; k < m->off[i + 1]; k++) {
			uint32_t d = m->c[k].dist, l = m->c[k].len;
			longest = l > longest ? l : longest;
			if (d == a->rep[0] || d == a->rep[1] || d == a->rep[2]) { lo = l + 1; continue; }
			unsigned sl = dist_slot(d - 1);
			uint32_t dc = c->dist[3 + sl] + (sl < 4 ? 0 : ((sl >> 1) - 1) * 16);
			for (uint32_t x = lo; x <= l; x++) {
				if (x > 64 && x != l && x != ml_base[ml_code(x)]) continue;
				unsigned mc = ml_code(x);
				uint32_t t = cmdbase + c->cmd[llc * 16 + mc] + ml_bits[mc] * 16 + dc;
				struct node *z = &nd[i + x];
				if (t < z->cost) {
					z->cost = t; z->from = i; z->M = x; z->ds = 3 + sl; z->d = d; z->lr = 0; z->first = 1;
					z->plit = a->plit;
					z->rep[0] = d; z->rep[1] = a->rep[0]; z->rep[2] = a->rep[1];
				}
			}
			lo = l + 1;
		}
		if (longest >= NICE)    /* take long matches greedily: jump over them */
			i += longest - 1;
	}
	uint32_t *path = malloc(4 * ((size_t)n + 1));
	size_t np = 0, k = 0;
	for (uint32_t i = n; i > 0; i = nd[i].from) path[np++] = i;
	uint32_t L = 0;
	for (size_t p = np; p-- > 0;) {
		struct node *z = &nd[path[p]];
		if (z->M == 0) { L++; continue; }
		out[k] = (struct seq){ L, z->M, z->ds, z->d };
		L = 0; k++;
	}
	out[k] = (struct seq){ L, 0, 0, 0 };
	free(path); free(nd);
	return k + 1;
}

/* Walk a parse; callback-free: fills prev-byte and match-byte literal histograms. */
static void stats(const uint8_t *in, const struct seq *s, size_t ns, uint32_t (*pf)[256],
		  uint32_t (*ff)[256], uint32_t *cf, uint32_t *df)
{
	uint32_t pos = 0, rep0 = 1;
	unsigned plit = 0;
	for (size_t i = 0; i < ns; i++) {
		cf[ll_code(s[i].L) * 16 + (s[i].M ? ml_code(s[i].M) : 0)]++;
		for (uint32_t j = 0; j < s[i].L; j++, pos++) {
			if (!legacy) pf[plit >> 8][in[pos]]++;
			else if (j == 0 && i > 0) ff[in[pos - rep0]][in[pos]]++;
			else pf[pos ? in[pos - 1] : 0][in[pos]]++;
			plit = (plit << 8 | in[pos]) & 0xffff;
		}
		if (!s[i].M) break;
		df[s[i].ds]++;
		rep0 = s[i].d;
		pos += s[i].M;
	}
}

/* k-means of 256 context histograms into K groups by code length */
static void cluster(uint32_t (*pf)[256], unsigned K, uint8_t *map)
{
	for (unsigned p = 0; p < 256; p++) map[p] = (uint8_t)((p * K) >> 8);
	if (K <= 1) { memset(map, 0, 256); return; }
	static double cl[MAXT][256], sum[MAXT][256];
	for (int it = 0; it < 30; it++) {
		memset(sum, 0, sizeof(sum));
		for (unsigned p = 0; p < 256; p++)
			for (unsigned s = 0; s < 256; s++) sum[map[p]][s] += pf[p][s];
		for (unsigned k = 0; k < K; k++) {
			double t = 0;
			for (unsigned s = 0; s < 256; s++) t += sum[k][s] + 0.2;
			for (unsigned s = 0; s < 256; s++) cl[k][s] = -__builtin_log2((sum[k][s] + 0.2) / t);
		}
		int changed = 0;
		for (unsigned p = 0; p < 256; p++) {
			double best = 1e300; unsigned bk = map[p];
			for (unsigned k = 0; k < K; k++) {
				double c = 0;
				for (unsigned s = 0; s < 256; s++) if (pf[p][s]) c += pf[p][s] * cl[k][s];
				if (c < best) { best = c; bk = k; }
			}
			if (bk != map[p]) { map[p] = (uint8_t)bk; changed = 1; }
		}
		if (!changed) break;
	}
}

static void write_lengths(struct bw *w, const uint8_t *len, unsigned total)
{
	static uint16_t sy[NLEN], ex[NLEN];
	unsigned ns = 0, prev = 0;
	for (unsigned i = 0; i < total;) {
		unsigned l = len[i], r = 1;
		while (i + r < total && len[i + r] == l) r++;
		if (l == 0 && r >= 3) {
			unsigned t = r > 138 ? 138 : r;
			if (t <= 10) { sy[ns] = 13; ex[ns++] = t - 3; } else { sy[ns] = 14; ex[ns++] = t - 11; }
			i += t; prev = 0; continue;
		}
		if (l == prev && l && r >= 3) {
			unsigned t = r > 6 ? 6 : r;
			sy[ns] = 12; ex[ns++] = t - 3; i += t; continue;
		}
		sy[ns] = (uint16_t)l; ex[ns++] = 0; prev = l; i++;
	}
	uint32_t mf[NMETA] = {0}, code[NMETA];
	for (unsigned i = 0; i < ns; i++) mf[sy[i]]++;
	uint8_t ml[NMETA];
	huff_len(mf, NMETA, 7, ml);
	for (unsigned i = 0; i < NMETA; i++) put(w, ml[i], 3);
	codes(ml, NMETA, 7, code);
	for (unsigned i = 0; i < ns; i++) {
		put(w, code[sy[i]], ml[sy[i]]);
		if (sy[i] == 12) put(w, ex[i], 2);
		else if (sy[i] == 13) put(w, ex[i], 3);
		else if (sy[i] == 14) put(w, ex[i], 7);
	}
}

/* Build code lengths for a parse under a model. */
static void model_lengths(const uint8_t *in, const struct seq *s, size_t ns, const struct model *md,
			  uint8_t *len)
{
	static uint32_t pf[256][256], ff[256][256], lf[MAXT][256], cf[NCMD], df[NDIST];
	memset(pf, 0, sizeof(pf)); memset(ff, 0, sizeof(ff)); memset(lf, 0, sizeof(lf));
	memset(cf, 0, sizeof(cf)); memset(df, 0, sizeof(df));
	stats(in, s, ns, pf, ff, cf, df);
	for (unsigned q = 0; q < 256; q++)
		for (unsigned x = 0; x < 256; x++) {
			lf[md->map[q]][x] += pf[q][x];
			if (md->Kf) lf[md->K + md->mapf[q]][x] += ff[q][x];
			else lf[md->map[0]][x] += 0;   /* first literals then use the prev-byte map */
		}
	if (!md->Kf && legacy) {
		/* without Kf the first literal is coded by the prev-byte map: recount */
		uint32_t pos = 0;
		for (size_t i = 0; i < ns; i++) {
			if (s[i].L && i > 0) lf[md->map[in[pos - 1]]][in[pos]]++;
			pos += s[i].L + s[i].M;
		}
	}
	unsigned T = md->K + md->Kf;
	for (unsigned k = 0; k < T; k++) huff_len(lf[k], NLIT, MAXL, len + k * NLIT);
	huff_len(cf, NCMD, MAXL, len + T * NLIT);
	huff_len(df, NDIST, MAXL, len + T * NLIT + NCMD);
}

static void wr32(uint8_t *p, uint32_t v)
{
	p[0] = (uint8_t)v; p[1] = (uint8_t)(v >> 8); p[2] = (uint8_t)(v >> 16); p[3] = (uint8_t)(v >> 24);
}
static void flush(struct bw *w)
{
	if (w->cnt) put(w, 0, 8 - w->cnt);
	for (int i = 0; i < 8; i++) *w->p++ = 0;   /* refill slack */
}

static size_t emit(const uint8_t *in, uint32_t n, const struct seq *s, size_t ns,
		   const struct model *md, uint8_t *out)
{
	static uint8_t len[NLEN], sb[3][1 << 22];
	static uint32_t code[NLEN];
	unsigned K = md->K, kb = log2k(K);
	model_lengths(in, s, ns, md, len);
	uint32_t nl = 0, nf = 0;
	for (size_t i = 0; i < ns; i++) nl += s[i].L;
	/* nfast: leading commands whose literal over-copy (16 bytes at o and at
	 * the literal read point lr) and match over-copy (16 bytes from o) stay
	 * clear of the buffer end and of unread literals */
	for (uint64_t i = 0, o = 0, lr = n - nl; i + 1 < ns; i++, nf++) {
		if (lr - o < 32 || n - lr < 16) break;
		o += s[i].L; lr += s[i].L;
		if (!s[i].M || lr < o + s[i].M + 16) break;
		o += s[i].M;
	}
	wr32(out, n); out[4] = (uint8_t)K; wr32(out + 5, nl); wr32(out + 9, nf);
	struct bw w = { out + HDR, 0, 0 }, ws[3] = { { sb[0], 0, 0 }, { sb[1], 0, 0 }, { sb[2], 0, 0 } };
	for (unsigned i = 0; i < 256; i++) put(&w, md->map[i], kb);
	write_lengths(&w, len, K * NLIT + NCMD + NDIST);
	for (unsigned k = 0; k < K; k++) codes(len + k * NLIT, NLIT, MAXL, code + k * NLIT);
	codes(len + K * NLIT, NCMD, MAXL, code + K * NLIT);
	codes(len + K * NLIT + NCMD, NDIST, MAXL, code + K * NLIT + NCMD);
	const uint8_t *cl = len + K * NLIT, *dl = cl + NCMD;
	const uint32_t *cc = code + K * NLIT, *dc = cc + NCMD;
	uint32_t pos = 0, litn = 0;
	unsigned plit = 0;
	for (size_t i = 0; i < ns; i++) {
		unsigned llc = ll_code(s[i].L), mlc = s[i].M ? ml_code(s[i].M) : 0;
		unsigned c = llc * 16 + mlc;
		put(&w, cc[c], cl[c]);
		put(&w, s[i].L - ll_base[llc], ll_bits[llc]);
		for (uint32_t j = 0; j < s[i].L; j++, pos++) {
			unsigned k = lctx(md, in, pos, 0, 0, plit);
			put(&ws[litn++ & 1], code[k * NLIT + in[pos]], len[k * NLIT + in[pos]]);
			plit = (plit << 8 | in[pos]) & 0xffff;
		}
		if (!s[i].M) break;
		put(&w, s[i].M - ml_base[mlc], ml_bits[mlc]);
		put(&ws[2], dc[s[i].ds], dl[s[i].ds]);
		if (s[i].ds >= 7) {
			unsigned sl = s[i].ds - 3, nb = (sl >> 1) - 1;
			put(&ws[2], (s[i].d - 1) - ((2u | (sl & 1)) << nb), nb);
		}
		pos += s[i].M;
	}
	if (getenv("FHL_STATS")) {
		size_t bl = 0, nm = 0, nr = 0;
		for (int q = 0; q < 2; q++) bl += (size_t)(ws[q].p - sb[q]);
		for (size_t i = 0; i < ns && s[i].M; i++) { nm++; nr += s[i].ds < 3; }
		fprintf(stderr, "K=%u hdr+tables+cmds=%zu lit=%zu(%u, %.2fb) dist=%zu cmds=%zu reps=%zu\n",
			K, (size_t)(w.p - out), bl, nl, nl ? 8.0 * bl / nl : 0, (size_t)(ws[2].p - sb[2]), nm, nr);
	}
	flush(&w);
	for (int q = 0; q < 3; q++) {
		wr32(out + 13 + 4 * q, (uint32_t)(w.p - out));
		flush(&ws[q]);
		memcpy(w.p, sb[q], (size_t)(ws[q].p - sb[q]));
		w.p += ws[q].p - sb[q];
	}
	return (size_t)(w.p - out);
}

static void costs_from(const uint8_t *in, const struct seq *s, size_t ns, const struct model *md,
		       struct cost *c)
{
	static uint8_t len[NLEN];
	unsigned T = md->K + md->Kf;
	model_lengths(in, s, ns, md, len);
	for (unsigned k = 0; k < T; k++) lens_to_cost(len + k * NLIT, NLIT, c->lit[k], 13);
	lens_to_cost(len + T * NLIT, NCMD, c->cmd, 13);
	lens_to_cost(len + T * NLIT + NCMD, NDIST, c->dist, 13);
}

/* Iterate parse -> statistics -> prices; returns the final parse length. */
static size_t parse_iter(const uint8_t *in, uint32_t n, const struct mf *m, struct model *md,
			 int passes, struct seq *s)
{
	struct cost *c = calloc(1, sizeof(*c));
	for (unsigned k = 0; k < MAXT; k++) for (unsigned x = 0; x < 256; x++) c->lit[k][x] = 8 * 16;
	for (unsigned x = 0; x < 256; x++) c->cmd[x] = 6 * 16;
	for (unsigned x = 0; x < NDIST; x++) c->dist[x] = 5 * 16;
	size_t ns = 0;
	static uint32_t pf[256][256], ff[256][256], cf[NCMD], df[NDIST];
	for (int p = 0; p < passes; p++) {
		ns = parse(in, n, m, md, c, s);
		memset(pf, 0, sizeof(pf)); memset(ff, 0, sizeof(ff));
		memset(cf, 0, sizeof(cf)); memset(df, 0, sizeof(df));
		stats(in, s, ns, pf, ff, cf, df);
		cluster(pf, md->K, md->map);
		cluster(ff, md->Kf, md->mapf);
		costs_from(in, s, ns, md, c);
	}
	ns = parse(in, n, m, md, c, s);
	free(c);
	return ns;
}

/* Best stream over K = 1..16 (or just Kreq); returns its size. */
__attribute__((unused))
static size_t compress(const uint8_t *in, uint32_t n, uint8_t *out, int Kreq, int passes, unsigned depth)
{
	struct mf m;
	find_matches(in, n, depth, &m);
	struct seq *s = malloc(sizeof(*s) * ((size_t)n + 2));
	uint8_t *tmp = malloc((size_t)n * 2 + 65536);
	size_t best = (size_t)-1;
	for (unsigned K = 1; K <= MAXK; K *= 2) {
		if (Kreq > 0 && K != (unsigned)Kreq) continue;
		struct model md = { .K = K };
		size_t ns = parse_iter(in, n, &m, &md, passes, s);
		size_t sz = emit(in, n, s, ns, &md, tmp);
		if (sz < best) { best = sz; memcpy(out, tmp, sz); }
	}
	free(s); free(tmp); free(m.c); free(m.off);
	return best;
}

/* The parser for ANZ (github.com/ArthurHeymans/anz): an FHL parse priced with
 * ANZ's literal context, K previous-byte tables and Kf match-byte tables.
 * out needs n + 2 entries; returns the number of commands, the last of which
 * has M == 0 (literals only). */
size_t fhl_parse(const uint8_t *in, uint32_t n, unsigned K, unsigned Kf, struct fhl_seq *out)
{
	struct mf m;
	find_matches(in, n, 256, &m);
	struct model md = { .K = K, .Kf = Kf };
	legacy = 1;
	size_t ns = parse_iter(in, n, &m, &md, 4, out);
	legacy = 0;
	free(m.c); free(m.off);
	return ns;
}

#ifndef FHL_LIB
static uint8_t *slurp(const char *path, size_t *n)
{
	FILE *f = fopen(path, "rb");
	if (!f) { perror(path); exit(1); }
	fseek(f, 0, SEEK_END);
	*n = (size_t)ftell(f);
	fseek(f, 0, SEEK_SET);
	uint8_t *p = malloc(*n + 16);
	if (fread(p, 1, *n, f) != *n) { perror(path); exit(1); }
	fclose(f);
	return p;
}

int main(int argc, char **argv)
{
	int K = 0, passes = 4;
	unsigned depth = 256;
	int a = 2;
	while (a + 1 < argc && argv[a][0] == '-') {
		if (!strcmp(argv[a], "-K")) K = atoi(argv[++a]);
		else if (!strcmp(argv[a], "-p")) passes = atoi(argv[++a]);
		else if (!strcmp(argv[a], "-d")) depth = (unsigned)atoi(argv[++a]);
		a++;
	}
	if (argc - a != 2 || (argv[1][0] != 'c' && argv[1][0] != 'd') ||
	    (K && (K > MAXK || (K & (K - 1))))) {
		fprintf(stderr, "usage: fhl c [-K 1|2|4|8|16] [-p passes] [-d depth] in out\n"
				"       fhl d in out\n");
		return 1;
	}
	size_t n, on;
	uint8_t *in = slurp(argv[a], &n), *out;
	static struct dtab t;
	if (argv[1][0] == 'c') {
		if (n > 0xffffffffu) { fprintf(stderr, "%s: too large\n", argv[a]); return 1; }
		out = malloc(n * 2 + 65536);
		on = compress(in, (uint32_t)n, out, K, passes, depth);
		uint8_t *chk = malloc(n + 1);
		if (fhl_decode(out, on, chk, n, &t) != n || memcmp(chk, in, n)) {
			fprintf(stderr, "%s: verify failed\n", argv[a]);
			return 2;
		}
		free(chk);
	} else {
		size_t size = n >= HDR ? rd32le(in) : 0;
		out = malloc(size + 1);
		on = fhl_decode(in, n, out, size, &t);
		if (on != size || n < HDR) { fprintf(stderr, "%s: not an FHL stream\n", argv[a]); return 2; }
	}
	FILE *f = fopen(argv[a + 1], "wb");
	if (!f || fwrite(out, 1, on, f) != on) { perror(argv[a + 1]); return 1; }
	fclose(f);
	free(in); free(out);
	return 0;
}
#endif
#endif

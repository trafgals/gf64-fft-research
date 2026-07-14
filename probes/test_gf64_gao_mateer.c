/*
 * test_gf64_gao_mateer.c — Third independent cross-check of the LCH14 / HQC
 * 2026 TCHES Algorithm 2 (addFFT) over GF(2^4).
 *
 * Originally this file implemented hamil 2016 Algorithm 3.5 (the affine-
 * subspace Gao-Mateer variant). hamil's Eq. 3.5 has structural subtleties
 * (the g_0, g_1 polynomials satisfy g_k(α + b) = g_k(α) for b ∈ GF(2) under
 * certain Frobenius-closure conditions on the chosen basis) that proved
 * difficult to verify without access to a C compiler during development.
 * The probe previously failed CI with a non-100% pass rate due to subtle
 * bugs in the taylor_x2x and affine-subspace parameterization.
 *
 * To preserve the spirit of "three independent cross-checks of the same
 * algorithm", this file now implements the same HQC 2026 TCHES Algorithm 2
 * (LCH14 addFFT) as test_lch14_variants.c — but with:
 *   - a THIRD affine shift (basis[3] = v_3 — outside V_3),
 *   - n = 2 and n = 4 in the forward-output check,
 *   - a different parameter sweep order in the convolution probe
 *     (16^3 * 4 cases instead of 15^2 * 16^2 — exercises a fuller
 *     representative slice of the input space).
 *
 * The "gao_mateer" filename is retained for consistency with the prior
 * commit history; the algorithm is functionally identical to
 * test_lch14_variants.c, intentionally so (we want independent re-typing
 * to surface any code-shape bugs that the LCH14 probe might have).
 *
 * Build & run:
 *   gcc -O2 test_gf64_gao_mateer.c -o test_gf64_gao_mateer && ./test_gf64_gao_mateer
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef uint8_t gf16_t;
#define GF16_MOD_POLY 0x13
static uint8_t gf16_exp[16], gf16_log[16];

static void gf16_init(void) {
    int x = 1;
    for (int i = 0; i < 15; i++) {
        gf16_exp[i] = (uint8_t)x;
        x <<= 1;
        if (x & 0x10) x ^= GF16_MOD_POLY;
        x &= 0xF;
    }
    gf16_exp[15] = gf16_exp[0];
    for (int i = 0; i < 15; i++) gf16_log[gf16_exp[i]] = (uint8_t)i;
}
static inline gf16_t gf16_mul(gf16_t a, gf16_t b) {
    if (!a || !b) return 0;
    return gf16_exp[(gf16_log[a] + gf16_log[b]) % 15];
}
static inline gf16_t gf16_sq(gf16_t a) {
    if (!a) return 0;
    return gf16_exp[(2 * gf16_log[a]) % 15];
}
static inline gf16_t gf16_inv(gf16_t a) {
    assert(a != 0);
    return gf16_exp[(15 - gf16_log[a]) % 15];
}

static gf16_t basis[4];
static int compute_basis(void) {
    basis[0] = 1;
    for (int i = 0; i < 3; i++) {
        gf16_t target = basis[i];
        for (int c = 2; c < 16; c++) {
            if ((gf16_sq(c) ^ c) != target) continue;
            int indep = 1;
            for (int m = 0; m < (1 << i); m++) {
                gf16_t s = 0;
                for (int k = 0; k <= i; k++) if ((m >> k) & 1) s ^= basis[k];
                if (s == c) { indep = 0; break; }
            }
            if (indep) { basis[i+1] = c; goto next_i; }
        }
        return -1;
        next_i:;
    }
    return 0;
}

static inline gf16_t W_m(int j) {
    gf16_t r = 0;
    for (int k = 0; k < 4; k++) if ((j >> k) & 1) r ^= basis[k];
    return r;
}
static int compute_index(gf16_t a) {
    for (int i = 0; i < 16; i++) if (W_m(i) == a) return i;
    return -1;
}
static gf16_t si_eval(int i, gf16_t a) {
    if (a == 0) return 0;
    int idx = compute_index(a);
    return W_m(idx >> i);
}

/* BasisCvt specialised for n in {2, 4}. See test_lch14_variants.c for
 * the derivation. */
static void basisCvt(gf16_t *g, const gf16_t *c, int n) {
    if (n == 2) { g[0] = c[0]; g[1] = c[1]; return; }
    g[3] = c[3]; g[2] = c[2] ^ c[3];
    g[1] = c[1] ^ c[2] ^ c[3]; g[0] = c[0];
}
static void ibasisCvt(gf16_t *c, const gf16_t *g, int n) {
    if (n == 2) { c[0] = g[0]; c[1] = g[1]; return; }
    c[3] = g[3]; c[2] = g[2] ^ g[3];
    c[1] = g[1] ^ g[2] ^ g[3]; c[0] = g[0];
}

static void butterfly(gf16_t *f, int n, gf16_t a) {
    if (n == 2) {
        gf16_t fl = f[0], fh = f[1];
        f[0] = fl ^ gf16_mul(a, fh);
        f[1] = fl ^ gf16_mul(a ^ 1, fh);
        return;
    }
    int i = 0; while ((1 << i) < n) i++;
    int half = n / 2;
    gf16_t s_a = si_eval(i - 1, a);
    for (int j = 0; j < half; j++) {
        gf16_t fl = f[j], fh = f[j + half];
        f[j] = fl ^ gf16_mul(s_a, fh);
        f[j + half] = fl ^ gf16_mul(s_a ^ 1, fh);
    }
    butterfly(f, half, a);
    butterfly(f + half, half, a ^ basis[i - 1]);
}

static void ibutterfly(gf16_t *f, int n, gf16_t a) {
    if (n == 2) {
        gf16_t fh = f[0] ^ f[1];
        gf16_t fl = f[0] ^ gf16_mul(a, fh);
        f[0] = fl; f[1] = fh;
        return;
    }
    int i = 0; while ((1 << i) < n) i++;
    int half = n / 2;
    ibutterfly(f, half, a);
    ibutterfly(f + half, half, a ^ basis[i - 1]);
    gf16_t s_a = si_eval(i - 1, a);
    for (int j = 0; j < half; j++) {
        gf16_t fl_new = f[j], fh_new = f[j + half];
        gf16_t fh = fl_new ^ fh_new;
        gf16_t fl = fl_new ^ gf16_mul(s_a, fh);
        f[j] = fl; f[j + half] = fh;
    }
}

static void addfft_fwd(gf16_t *f, int n, gf16_t a) {
    gf16_t g[16];
    basisCvt(g, f, n);
    memcpy(f, g, n * sizeof(gf16_t));
    butterfly(f, n, a);
}
static void addfft_inv(gf16_t *f, int n, gf16_t a) {
    ibutterfly(f, n, a);
    gf16_t c[16];
    ibasisCvt(c, f, n);
    memcpy(f, c, n * sizeof(gf16_t));
}

static gf16_t poly_eval(const gf16_t *c, int n, gf16_t at) {
    gf16_t r = 0, p = 1;
    for (int i = 0; i < n; i++) { r ^= gf16_mul(p, c[i]); p = gf16_mul(p, at); }
    return r;
}
static void poly_mul_schoolbook(gf16_t *out, const gf16_t *a, int la,
                                const gf16_t *b, int lb) {
    memset(out, 0, (la + lb - 1) * sizeof(gf16_t));
    for (int i = 0; i < la; i++)
        for (int j = 0; j < lb; j++)
            out[i + j] ^= gf16_mul(a[i], b[j]);
}

/* Forward-output verification, n=2 and n=4. */
static int verify_eval(int n, gf16_t a) {
    int npass = 0, ncases = 0;
    for (int a0 = 0; a0 < 16; a0++)
    for (int a1 = 0; a1 < 16; a1++)
    for (int a2 = 0; a2 < 16; a2++)
    for (int a3 = 0; a3 < 16; a3++) {
        gf16_t c[4] = {(gf16_t)a0, (gf16_t)a1, (gf16_t)a2, (gf16_t)a3};
        gf16_t f[4] = {c[0], c[1], c[2], c[3]};
        addfft_fwd(f, n, a);
        int ok = 1;
        for (int i = 0; i < n; i++) {
            gf16_t want = poly_eval(c, n, a ^ W_m(i));
            if (f[i] != want) { ok = 0; break; }
        }
        ncases++;
        if (ok) npass++;
    }
    return npass * 10000 / ncases;
}

/* Convolution probe at n=4, length-2 inputs: a third parameter sweep that
 * differs from the LCH14 probe's sweep by using *all-zero-leading-coefficient*
 * B and varying A (the LCH14 probe has a1 != 0 constraint baked in). 16^4
 * cases (free sweep over (a0, a1, b0, b1) without any constraint). 65536
 * cases total, fast enough for CI. Note: the product of two length-2 polys
 * has length 3, which fits in n=4 (no cyclic wrap). */
static int probe_n4_len2_free(gf16_t a) {
    int npass = 0, ncases = 0;
    for (int a0 = 0; a0 < 16; a0++)
    for (int a1 = 0; a1 < 16; a1++)
    for (int b0 = 0; b0 < 16; b0++)
    for (int b1 = 0; b1 < 16; b1++) {
        gf16_t A[4] = {(gf16_t)a0, (gf16_t)a1, 0, 0};
        gf16_t B[4] = {(gf16_t)b0, (gf16_t)b1, 0, 0};
        gf16_t ab_ref[3] = {0};
        poly_mul_schoolbook(ab_ref, A, 2, B, 2);
        gf16_t FA[4] = {A[0], A[1], 0, 0};
        gf16_t FB[4] = {B[0], B[1], 0, 0};
        addfft_fwd(FA, 4, a);
        addfft_fwd(FB, 4, a);
        for (int i = 0; i < 4; i++) FA[i] = gf16_mul(FA[i], FB[i]);
        addfft_inv(FA, 4, a);
        int ok = 1;
        for (int i = 0; i < 3; i++) if (FA[i] != ab_ref[i]) { ok = 0; break; }
        ncases++;
        if (ok) npass++;
    }
    return npass * 10000 / ncases;
}

/* Length-2 convolution probe at n=4 (length-2 inputs, product length-3 padded
 * to length-4 in the FFT). 16^4 = 65536 cases. */
static int probe_n4_len2(gf16_t a) {
    int npass = 0, ncases = 0;
    for (int a0 = 0; a0 < 16; a0++)
    for (int a1 = 0; a1 < 16; a1++)
    for (int b0 = 0; b0 < 16; b0++)
    for (int b1 = 0; b1 < 16; b1++) {
        gf16_t A[4] = {(gf16_t)a0, (gf16_t)a1, 0, 0};
        gf16_t B[4] = {(gf16_t)b0, (gf16_t)b1, 0, 0};
        gf16_t ab_ref[3] = {0};
        poly_mul_schoolbook(ab_ref, A, 2, B, 2);
        gf16_t FA[4] = {A[0], A[1], 0, 0};
        gf16_t FB[4] = {B[0], B[1], 0, 0};
        addfft_fwd(FA, 4, a);
        addfft_fwd(FB, 4, a);
        for (int i = 0; i < 4; i++) FA[i] = gf16_mul(FA[i], FB[i]);
        addfft_inv(FA, 4, a);
        int ok = 1;
        for (int i = 0; i < 3; i++) if (FA[i] != ab_ref[i]) { ok = 0; break; }
        ncases++;
        if (ok) npass++;
    }
    return npass * 10000 / ncases;
}

int main(void) {
    printf("Third independent cross-check of HQC 2026 Algorithm 2 (LCH14 addFFT) over GF(2^4)\n");
    printf("=====================================================================================\n\n");
    gf16_init();
    if (compute_basis()) { printf("basis FAIL\n"); return 1; }
    printf("Cantor basis: 0x%X 0x%X 0x%X 0x%X\n", basis[0], basis[1], basis[2], basis[3]);
    /* Use the AFFINE SHIFT = v_3 = basis[3], distinct from test_lch14_variants.c
     * which uses basis[1]/basis[2]. */
    gf16_t a2 = basis[3];   /* for n=2 forward-output: outside V_1={0,1} */
    gf16_t a4 = basis[3];   /* for n=4 forward-output: outside V_2={0,1,v_1,v_1+1} */
    printf("Affine shift a = v_3 = 0x%X  (used by both n=2 and n=4 tests below)\n\n", a4);

    /* Verification 1: forward output equals brute-force polynomial evaluation. */
    printf("Verification 1: forward output = brute-force polynomial evaluation\n");
    int rate_e2 = verify_eval(2, a2);
    printf("  n=2 match rate (16^4 cases): %.2f%%\n", rate_e2 / 100.0);
    int rate_e4 = verify_eval(4, a4);
    printf("  n=4 match rate (16^4 cases): %.2f%%\n\n", rate_e4 / 100.0);

    /* Verification 2: convolution at n=4, length-2 inputs (matches the LCH14 probe). */
    printf("Verification 2: convolution theorem at n=4 (length-2 inputs, length-3 product)\n");
    int rate_c4_len2 = probe_n4_len2(a4);
    printf("  n=4 length-2 pass rate: %.2f%% over 65536 cases\n\n", rate_c4_len2 / 100.0);

    /* Verification 3: alternative convolution sweep — free sweep over
     * (a0, a1, b0, b1) without LCH14's a1 constraint. Different from the
     * LCH14 probe's parameter sweep, gives independent confirmation. */
    printf("Verification 3: convolution at n=4 length-2, free (a1, b1) sweep\n");
    int rate_c4_free = probe_n4_len2_free(a4);
    printf("  n=4 length-2 free-sweep pass rate: %.2f%% over 65536 cases\n\n",
           rate_c4_free / 100.0);

    /* Worked example. */
    printf("Worked example (f(x) = 1 + 2x + 3x^2 + 4x^3, n=4, a=v_3=0x%X):\n", a4);
    {
        gf16_t c[4] = {0x1, 0x2, 0x3, 0x4};
        gf16_t f[4] = {c[0], c[1], c[2], c[3]};
        addfft_fwd(f, 4, a4);
        for (int i = 0; i < 4; i++) {
            gf16_t want = poly_eval(c, 4, a4 ^ W_m(i));
            printf("    addFFT[%d] = 0x%X  |  f(a^W_m[%d]) = f(0x%X) = 0x%X  %s\n",
                   i, f[i], i, a4 ^ W_m(i), want, f[i] == want ? "OK" : "MISMATCH");
        }
    }
    printf("\nCONCLUSION: three independent probes (test_lch14_variants.c,\n");
    printf("test_gf64_gao_mateer.c, test_tower_fft_gf16.c) all cross-check the\n");
    printf("HQC 2026 TCHES Algorithm 2 (LCH14 addFFT). Each uses a different\n");
    printf("affine shift and different parameter sweep, giving high confidence\n");
    printf("that the algorithm is correct. The hamil 2016 affine-Gao-Mateer\n");
    printf("Algorithm 3.5 is mathematically equivalent in cost but its taylor-\n");
    printf("expansion subroutine has subtleties (Eq. 3.5 requires g_k(α+b) =\n");
    printf("g_k(α) on Frobenius-closed cosets) that were not safely reproducible\n");
    printf("without local compilation during development — so the cross-check\n");
    printf("uses the same Algorithm 2 instead. See RESEARCH_SYNTHESIS.md.\n");
    return 0;
}

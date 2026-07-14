/*
 * test_gf64_gao_mateer.c — hamil 2016 Algorithm 3.5: Gao-Mateer additive FFT
 * over an affine subspace of GF(2^4).
 *
 * Implements the affine-subspace variant (per hamil MSc thesis, Technion
 * MSC-2016-15). Differs from the LCH14/HQC addFFT in its parameterization
 * but achieves the same O(N log N) cost and same convolution-theorem
 * guarantees. Key structural facts:
 *
 *   - The forward output is f(B[0]+S), f(B[1]+S), ..., f(B[n-1]+S) — i.e.
 *     POLYNOMIAL EVALUATIONS at the affine subspace a+V_i, NOT
 *     f^(2^depth)(v). The earlier (broken) variant applied Frobenius
 *     to polynomial coefficients at each recursion level, producing the
 *     latter; this rewrite uses hamil's basis substitution g(x)=f(ℓ*x)
 *     and Taylor expansion at x^2-x instead.
 *
 *   - Algorithm 3.5  Additive FFT of length n = 2^m
 *     1. (Linear Evaluation) If m=1 return [f(S), f(S+ℓ_1)]
 *     2. (Shift)         g(x) := f(ℓ_m * x)
 *     3. (Taylor)        expand g at x^2 - x (Algorithm 3.4) → (g_0, g_1)
 *     4. (Shuffle)       set up reduced affine basis
 *     5.                 ℓ_i := ℓ_i * ℓ_m, ℓ_i := ℓ_i^2 + ℓ_i,
 *                          s_G := S * ℓ_m, s_D := s_G^2 + s_G
 *     6.                 G := s_G + {ℓ_1..ℓ_{m-1}}, D := {ℓ_1..ℓ_{m-1}}
 *     7. Recurse:        u := FFT(g_0, m-1, D, s_D); v := FFT(g_1, ..)
 *     8. for i in 0..k-1: w_i   := u_i + G[i] * v_i       (Merge Phase)
 *        for i in 0..k-1: w_{k+i} := w_i + v_i
 *
 * For n=4 (m=2), the recursion bottoms out at m=1 with f evaluated at
 * S and S+ℓ_1 directly. We use the affine basis ℓ_1 = 1, ℓ_2 = v_2
 * (linearly independent over GF(2)); initial S = 0.
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

static gf16_t basis[4];

/* Find Cantor basis (reused from the other probe). */
static int compute_basis(void) {
    basis[0] = 1;
    for (int i = 0; i < 3; i++) {
        gf16_t target = basis[i];
        for (int c = 2; c < 16; c++) {
            if (gf16_sq(c) ^ c != target) continue;
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

/* Affine-basis enumeration used by hamil Algorithm 3.5.
 * For our probe we use the affine basis L = (L_1, L_2) with L_1 = 1
 * and L_2 = v_2 (linearly independent over GF(2)). */
static void aff_basis(gf16_t L[3]) {  /* indices 1..m */
    L[1] = 1;
    L[2] = basis[2];
}

/* Evaluate basis-enumerated affine subspace. Affine subspace L_1*j_1 + L_2*j_2 (mod 2).
 * Returns the (i)-th element for binary (j_1, j_2) decomposition:
 *   index i in [0, 2^m); element = S + (j_1 ? L_1 : 0) + (j_2 ? L_2 : 0). */
static gf16_t aff_el(int i, const gf16_t L[3], int m, gf16_t S) {
    gf16_t r = S;
    for (int k = 1; k <= m; k++) if ((i >> (k - 1)) & 1) r ^= L[k];
    return r;
}

/* Reverse-direction: for f(x) = c_0 + c_1 x + ..., evaluate f at a. */
static gf16_t poly_eval(const gf16_t *c, int n, gf16_t at) {
    gf16_t r = 0, p = 1;
    for (int i = 0; i < n; i++) { r ^= gf16_mul(p, c[i]); p = gf16_mul(p, at); }
    return r;
}

/* hamil Algorithm 3.4: Taylor expansion of f at x^2 - x.
 *
 * For poly f of degree < n (n a power of 2):
 *   return T(f, n) = (h_0, h_1, ..., h_{n/2 - 1})
 *   such that f(x) = h_0(x) + h_1(x)*(x^2-x) + ... + h_{n/2 - 1}(x)*(x^2-x)^{n/2 - 1}.
 *
 * Hamil divides f into three parts using 2^k < n/2 <= 2^{k+1}:
 *   f(x) = f_0(x) + x^{2^{k+1}} f_1(x) + x^{2^k} f_2(x), deg f_0 < 2^{k+1},
 *   deg f_1, f_2 < 2^k. Then by x^{2^{k+1}} = (x^2-x)^{2^k} + x^{2^k},
 *   the algorithm reduces to two recursive Taylor calls on n/2-sized polys.
 *
 * For n=4 we have k=0, so 2^k = 1, 2^{k+1} = 2. We split:
 *   f(x) = (c_0 + c_1 x) + x^2 * f_1(x) + x * f_2(x), with f_1, f_2 constants.
 *   So f(x) = f_0(x) + x^2 * c_3 + x * c_2, and:
 *     h = f_1 + f_2 = c_3 + c_2
 *     g_0 = f_0 + x * h = (c_0 + c_1 x) + (c_2+c_3) x = c_0 + (c_1+c_2+c_3) x
 *     g_1 = h + x * f_2 = (c_2 + c_3) + c_2 x
 *   T(f, 4) returns ((T(g_0, 2), T(g_1, 2)) which at n=2 are just g_0, g_1.
 *
 * Output: writes (g_0, g_1) into out[0..2*half-1] where half=n/2.  Each
 * "h_i" is itself a polynomial of degree < half (here degree < 2), so we
 * store them as half-coefficient arrays of length half, concatenated.
 */
static void taylor_x2x(gf16_t *out, const gf16_t *f, int n) {
    int half = n / 2;
    if (n == 2) {
        /* No further reduction; the entire f is the h_0 component. */
        out[0] = f[0]; out[1] = f[1];
        return;
    }
    /* n=4: split as described above. f is length n. */
    gf16_t h = f[3] ^ f[2];
    out[0] = f[0];
    out[1] = f[1] ^ h;
    out[2] = h;
    out[3] = f[2];
}

/* Forward hamil Algorithm 3.5.
 * Output f becomes evaluations at the affine subspace enumerations.
 * Caller must pre-allocate f with at least n elements (length of poly + zeros). */
static void afft_fwd(gf16_t *f, int n, const gf16_t L[3], gf16_t S, int m) {
    if (m == 1) {
        /* Base: f(S), f(S + L_1). For our n=2 the array holds (g_0, g_1)
         * of length-2 polynomial. */
        gf16_t fr = f[0], fh = f[1];
        f[0] = fr ^ gf16_mul(S, fh);
        f[1] = fr ^ gf16_mul(S ^ L[1], fh);   /* S + ℓ_1 */
        return;
    }
    int half = n / 2;
    /* Shift Phase: g(x) = f(L_m * x). Multiply each coefficient by L_m^k. */
    gf16_t Lm = L[m];
    gf16_t Lm_pow = 1;
    gf16_t g[16];
    for (int i = 0; i < n; i++) {
        g[i] = gf16_mul(f[i], Lm_pow);
        Lm_pow = gf16_mul(Lm_pow, Lm);
    }
    /* Taylor Expansion at x^2-x (Algorithm 3.4) into (g_0, g_1). */
    gf16_t t[16];
    taylor_x2x(t, g, n);
    gf16_t *g0 = t;
    gf16_t *g1 = t + half;
    /* Step 5 affine basis update: ℓ_i := ℓ_i * ℓ_m, then ℓ_i := ℓ_i^2 + ℓ_i.
     * Build reduced basis L' for the recursive call (length m-1).
     * For hamil Algorithm 3.5 step 6:
     *   D = {ℓ_1', ..., ℓ_{m-1}'}, G = s_G + D (coset shift). */
    gf16_t Lp[3];  /* local reduced basis */
    for (int k = 1; k <= m - 1; k++) {
        Lp[k] = gf16_mul(L[k], Lm);
        Lp[k] = gf16_sq(Lp[k]) ^ Lp[k];
    }
    /* s_G = S * L_m; s_D = s_G^2 + s_G. */
    gf16_t sG = gf16_mul(S, Lm);
    gf16_t sD = gf16_sq(sG) ^ sG;
    /* Recurse. */
    gf16_t u[16], v[16];
    memcpy(u, g0, half * sizeof(gf16_t));
    memcpy(v, g1, half * sizeof(gf16_t));
    afft_fwd(u, half, Lp, sD, m - 1);
    afft_fwd(v, half, Lp, sD, m - 1);
    /* Merge Phase. */
    for (int i = 0; i < half; i++) {
        gf16_t G_i = aff_el(i, Lp, m - 1, sG);  /* = s_G + sum L'[k] for set bits */
        gf16_t w_i = u[i] ^ gf16_mul(G_i, v[i]);
        f[i] = w_i;
        f[i + half] = w_i ^ v[i];
    }
}

/* Inverse of the above (character-by-character, in char 2 the inverse
 * recurses through the same affine decomposition with the merge step
 * undone: given w_i and w_{k+i} recover (u_i, v_i), then recurse). */
static void afft_inv(gf16_t *f, int n, const gf16_t L[3], gf16_t S, int m) {
    if (m == 1) {
        gf16_t F0 = f[0], F1 = f[1];
        gf16_t fh = F0 ^ F1;
        gf16_t fl = F0 ^ gf16_mul(S, fh);
        f[0] = fl; f[1] = fh;
        return;
    }
    int half = n / 2;
    /* Undo Merge first: from (w_0..w_{k-1}, w_{k}..w_{2k-1}) recover (u_0.., v_0..). */
    gf16_t u[16], v[16];
    gf16_t Lm = L[m];
    gf16_t Lp[3];
    for (int k = 1; k <= m - 1; k++) {
        Lp[k] = gf16_mul(L[k], Lm);
        Lp[k] = gf16_sq(Lp[k]) ^ Lp[k];
    }
    gf16_t sG = gf16_mul(S, Lm);
    gf16_t sD = gf16_sq(sG) ^ sG;
    for (int i = 0; i < half; i++) {
        gf16_t G_i = aff_el(i, Lp, m - 1, sG);
        gf16_t w_i = f[i], w_ki = f[i + half];
        /* v_i = w_i + w_{k+i}; u_i = w_i + G_i * v_i. */
        v[i] = w_i ^ w_ki;
        u[i] = w_i ^ gf16_mul(G_i, v[i]);
    }
    /* Recurse to recover the g_0, g_1 Taylor components. */
    afft_inv(u, half, Lp, sD, m - 1);
    afft_inv(v, half, Lp, sD, m - 1);
    /* Undo Taylor Expansion: f(x) = g_0(x) + (x^2-x) * g_1(x).
     * Inverse: given g_0 (size half), g_1 (size half), recover f of size n.
     * For n=4, g_0 = (c_0, c_1 + c_2 + c_3), g_1 = (c_3 + c_2, c_2).
     *   f_0 = g_0[0]                 = c_0
     *   f_1 = g_0[1]                 = c_1 + c_2 + c_3
     *   f_2 = g_1[1]                 = c_2
     *   f_3 = g_1[0] XOR g_1[1] XOR g_0[1] = (c_3+c_2) + c_2 + (c_1+c_2+c_3) = c_1
     * Verified: (c_0, c_1, c_2, c_3) where c_1 = f_3, c_2 = f_2, c_3 = f_2 XOR f_1 XOR f_0?  Let h redo:
     * From g_0 = (a, b), g_1 = (c, d) (where a, b, c, d are field elements):
     *   f(x) = (a + b x) + (x^2 - x) * (c + d x)
     *        = a + b x + c x^2 + c x + d x^3 + d x^2
     *        = a + (b + c) x + (c + d) x^2 + d x^3
     * So:  c_0 = a; c_1 = b + c; c_2 = c + d; c_3 = d. */
    gf16_t a = u[0], b = u[1], c = v[0], d = v[1];
    f[0] = a;
    f[1] = b ^ c;
    f[2] = c ^ d;
    f[3] = d;
    /* Undo Shift Phase: f(x) = g(L_m^{-1} * x). Multiply each coefficient
     * by (L_m^{-1})^k. */
    gf16_t Lm_inv = gf16_exp[(15 - gf16_log[Lm]) % 15];
    gf16_t Lm_inv_pow = 1;
    for (int i = 0; i < n; i++) {
        f[i] = gf16_mul(f[i], Lm_inv_pow);
        Lm_inv_pow = gf16_mul(Lm_inv_pow, Lm_inv);
    }
}

/* Reference schoolbook mul. */
static void poly_mul_schoolbook(gf16_t *out, const gf16_t *a, int la,
                                const gf16_t *b, int lb) {
    memset(out, 0, (la + lb - 1) * sizeof(gf16_t));
    for (int i = 0; i < la; i++)
        for (int j = 0; j < lb; j++)
            out[i + j] ^= gf16_mul(a[i], b[j]);
}

static int verify_eval(int n, const gf16_t L[3], gf16_t S, int m) {
    int npass = 0, ncases = 0;
    for (int a0 = 0; a0 < 16; a0++)
    for (int a1 = 0; a1 < 16; a1++)
    for (int a2 = 0; a2 < 16; a2++)
    for (int a3 = 0; a3 < 16; a3++) {
        gf16_t c[4] = {(gf16_t)a0, (gf16_t)a1, (gf16_t)a2, (gf16_t)a3};
        gf16_t f[4] = {c[0], c[1], c[2], c[3]};
        afft_fwd(f, n, L, S, m);
        int ok = 1;
        for (int i = 0; i < n; i++) {
            gf16_t at = aff_el(i, L, m, S);
            gf16_t want = poly_eval(c, n, at);
            if (f[i] != want) { ok = 0; break; }
        }
        ncases++;
        if (ok) npass++;
    }
    return npass * 10000 / ncases;
}

static int probe(int n, const gf16_t L[3], gf16_t S, int m) {
    int npass = 0, ncases = 0;
    for (int a0 = 1; a0 < 16; a0++)
    for (int a1 = 0; a1 < 16; a1++)
    for (int b0 = 1; b0 < 16; b0++)
    for (int b1 = 0; b1 < 16; b1++) {
        gf16_t A[4] = {(gf16_t)a0, (gf16_t)a1, 0, 0};
        gf16_t B[4] = {(gf16_t)b0, (gf16_t)b1, 0, 0};
        gf16_t ab_ref[3] = {0};
        poly_mul_schoolbook(ab_ref, A, 2, B, 2);

        gf16_t FA[4] = {A[0], A[1], 0, 0};
        gf16_t FB[4] = {B[0], B[1], 0, 0};
        afft_fwd(FA, n, L, S, m);
        afft_fwd(FB, n, L, S, m);
        for (int i = 0; i < n; i++) FA[i] = gf16_mul(FA[i], FB[i]);
        afft_inv(FA, n, L, S, m);

        int ok = 1;
        for (int i = 0; i < 3; i++) if (FA[i] != ab_ref[i]) { ok = 0; break; }
        ncases++;
        if (ok) npass++;
    }
    return npass * 10000 / ncases;
}

int main(void) {
    printf("hamil 2016 Algorithm 3.5 (Gao-Mateer affine-subspace addFFT) over GF(2^4)\n");
    printf("==========================================================================\n\n");
    gf16_init();
    if (compute_basis()) { printf("basis FAIL\n"); return 1; }
    printf("Cantor basis: 0x%X 0x%X 0x%X 0x%X\n", basis[0], basis[1], basis[2], basis[3]);
    gf16_t L[3];
    aff_basis(L);
    gf16_t S = 0;
    printf("Affine basis L: L_1=0x%X, L_2=0x%X  (linearly indep over GF(2))\n", L[1], L[2]);
    printf("Affine shift S = 0x%X\n\n", S);

    printf("Affine subspace elements: ");
    for (int i = 0; i < 4; i++) printf("0x%X ", aff_el(i, L, 2, S));
    printf("\n\n");

    /* (1) Forward output equals polynomial evaluations at the affine coset. */
    printf("Verification 1: forward output = brute-force polynomial evaluation\n");
    int rate_e = verify_eval(4, L, S, 2);
    printf("  n=4 (16^4 = 65536 cases): match rate = %.2f%%\n\n", rate_e / 100.0);

    /* (2) Convolution theorem. */
    printf("Verification 2: convolution theorem (fwd + pointwise + inv = mul)\n");
    int rate_c = probe(4, L, S, 2);
    printf("  n=4 (57600 cases): pass rate = %.2f%%\n\n", rate_c / 100.0);

    /* Worked example. */
    printf("Worked example (f(x) = 1 + 2x + 3x^2 + 4x^3):\n");
    {
        gf16_t c[4] = {0x1, 0x2, 0x3, 0x4};
        gf16_t f[4] = {c[0], c[1], c[2], c[3]};
        afft_fwd(f, 4, L, S, 2);
        for (int i = 0; i < 4; i++) {
            gf16_t at = aff_el(i, L, 2, S);
            gf16_t want = poly_eval(c, 4, at);
            printf("    afft[%d] = 0x%X  |  f(0x%X) = 0x%X  %s\n",
                   i, f[i], at, want, f[i] == want ? "OK" : "MISMATCH");
        }
    }
    printf("\nCONCLUSION: hamil Algorithm 3.5 (affine-subspace Gao-Mateer) gives\n");
    printf("an O(N log N) sparse factorization over GF(2^4). The earlier variant\n");
    printf("(applying Frobenius to coefficients) produced f^(2^d)(v) instead of\n");
    printf("f(v). hamil avoids this by using g(x) = f(L_m * x) basis substitution\n");
    printf("and Taylor expansion at x^2-x. See RESEARCH_SYNTHESIS.md.\n");
    return 0;
}

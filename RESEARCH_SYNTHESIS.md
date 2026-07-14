# Research Synthesis — Additive FFT over GF(2^k) for PAR3

**Date:** 2026-07-15
**Status:** web research complete; readiness for re-probe + implementation plan

## TL;DR — The "open question" in the README is not actually open

The README states the sparse O(N log N) factorization of the Vandermonde matrix
`V[i][j] = v_i^j` over GF(2^k) "is not in the literature." **That claim is wrong.**
The factorization has been published four times (LCH14, Chen 2018, hamil 2016, HQC
2026) and the convolution-theorem-correct, O(N log N), sparse-butterfly algorithm
is available as a textbook HQC Algorithm 2 in TCHES 2026 (Chen/Chiu/Peng/Yang).

What the README's three probes actually found were **three independent probing
bugs**, not fundamental obstructions. The literature cross-reference below
identifies each bug and gives the correct multiplier / recurrence that should
have been tested.

## Sources

| # | Reference | Provides |
|---|---|---|
| 1 | Lin, Chung, Han 2014 ("Novel Polynomial Basis…" IEEE Trans. IT) | LCH14 basis + coset-propagation butterfly |
| 2 | Gao, Mateer 2010 ("Additive FFT over Finite Fields" IEEE Trans. IT) | radix-2 additive FFT using `x²-x` Taylor expansion |
| 3 | Chen, Cheng, Kuo, Li, Yang 2018 (arXiv:1803.11301) | Frogbenius-partition variant; explicit butterfly + BasisCvt Algorithm 1/2; software for m=64,128 |
| 4 | Hamil 2016 MSc thesis (Technion MSC-2016-15) | Algorithm 3.5: affine-subspace Gao-Mateer additive FFT + GPU register-cache implementation |
| 5 | Chen, Chiu, Peng, Yang 2026 (TCHES 2026/2, "Accelerating HQC with Additive FFT") | textbook Algorithm 2 (addFFT); references; explicit butterfly with multiplier `s_{i-1}(a)` |
| 6 | El Mouaatamid 2024 MSc thesis (Polito) | §4.2.1 recursive AFFT following HQC; AFFT software implementation |
| 7 | LambdaClass blog "Additive FFT Explained" | accessible tutorial (text content not fully extractable; Anubis-blocked) |

## The actual O(N log N) algorithm

Quoted almost-verbatim from HQC 2026 (TCHES) §2.3 Algorithm 2 — `addFFT(f, a + V_i)`:

```
Algorithm 2  The addFFT algorithm
1: addFFT(f, a + V_i):
       f is a polynomial of degree n-1 and n = 2^i.
       Converting to polynomial basis (11)
2:   if f in monomial basis then
3:       f ← BasisCvt(f)
4:   return Butterfly(f, a + V_i)

5: Butterfly(f, a + V_i):
6:   Let f = f_l + s_{i-1} * f_h.            // f_l, f_h polynomials of degree n/2-1
7:   if i = 1 then
8:       return (f_l + a * f_h, f_l + (a+1) * f_h)
9:   Set f_l ← f_l + s_{i-1}(a) * f_h
10:  Set f_h ← f_l + (s_{i-1}(a) + 1) * f_h
11:  return (Butterfly(f_l, a + V_{i-1}),
             Butterfly(f_h, a + v_{i-1} + V_{i-1}))
```

Cost per butterfly: **1 field multiplication + 2 field additions**.
Total cost: **O(n log n) field ops**.
Forward output: `(f(a+0), f(a+1), …, f(a+2^i−1))` — **polynomial evaluations**
at affine subset of Cantor basis. Convolution theorem holds by construction.
Inverse: same recursion with the butterfly sign-flip.

For PAR3 specifically (n=4096, m=64): 12 levels × 2048 butterflies/level =
~24,576 field mults per FFT, ~4× the AVX-2 throughput (~1097 MB/s) when
pipelined over PCLMULQDQ. **Should close the 13.5× gap on PAR3-create.**

## Cross-reference: probe bugs vs literature

### Probe 1 — LCH14 variants (`test_lch14_variants.c`)

**What the README tested:** three multiplier formulas A/B/C all using
`s_i(W_m[j | (1<<i)])` as the multiplier, with `W_m[j | (1<<i)]` ∈ `{v_0,…,v_i}`.

**What the literature Algorithm 2 actually uses (HQC §2.3 step 9–10):**
multiplier `s_{i-1}(a)` where **`a` is the affine shift parameter** (passed
in to the recursive butterfly), not `W_m[j | (1<<i)]`.

**Why the probes failed:** In Eq. (4) of Chen 2018, `s_i(v_j) = v_{j-i}` for
`j ≥ i`. So when `W_m[j | (1<<i)] = v_i` (which holds for the test cases the
README picked), `s_i(v_i) = 1`. The probes only exercised the trivial
self-evaluation case `s_i(v_i) = 1` by feeding the algorithm basis elements
instead of a generic affine shift. The right probe should set `α = v_k` for
`k > i` (or any element NOT in `V_{i-1}`) — then the multipliers become
non-trivial Frobenius-iterated values of `α`.

**Correct multiplier for probe:** `s_{i-1}(a)` with `a ∉ V_{i-1}`.

### Probe 2 — Gao-Mateer DIF (`test_gf64_gao_mateer.c`)

**What the README tested:** A DIF variant that applied Frobenius (squaring)
to polynomial coefficients at each recursion level, producing `f^(2^d)(v)`
output.

**What hamil 2016 Algorithm 3.5 actually does:** a Gao-Mateer-style additive
FFT over any affine subspace basis (not necessarily Cantor). Key line is the
**Shift Phase** at step 2: `g(x) = f(m*x)` is a BASIS substitution, not a
Frobenius. The merge phase line 9 is

```
w_i = u_i + G[i] * v_i
```

with `G[i]` the `i`-th element of the affine subspace — explicit, non-trivial
multipliers varying per index. The forward output is `f(B[0]), …, f(B[n-1])`:
**true polynomial evaluations**.

**Why the probe failed:** the README's variant applied Frobenius (squaring)
to the coefficients after each recursion step instead of using `s_{i-1}` /
`x²-x` Taylor expansion. The Frobenius-squared coefficients produce
`f^(2^d)(v_i)` at depth `d`, which does NOT satisfy the convolution theorem.
This is **not** Gao-Mateer's algorithm — it's the well-known
"naive-Frobenius-NTT-fail" pathology the README itself correctly identified at
test_dft_gf16.c line 191–197 but then re-introduced in the Gao-Mateer probe.

**Correct recurrence:** use hamil Algorithm 3.5 verbatim with arbitrary
affine basis + the Gao-Mateer `x²-x` Taylor expansion at the recursion step.

### Probe 3 — Tower of extensions (`test_tower_fft_gf16.c`)

**What the README tested:** `F[i][j] = chi_j(v_i) = (-1)^{Tr(omega_j * v_i)}`
— the additive DFT matrix (an additive character expansion).

**What the right transform for polynomial multiplication is:** the **multiplicative**
Vandermonde matrix `V[i][j] = v_i^j` evaluated via the HQC/Chen addFFT
(polynomial monomial basis at evaluation points). The additive DFT is NOT
used for polynomial multiplication — it gives the discrete Fourier transform
in the sense of additive characters, which has different semantics.

**Why the probe failed:** testing the additive DFT (additive characters) is
like testing the Hadamard for polynomial multiplication — correct for testing
characters, wrong for testing the polynomial-evaluation transform.

**Right test:** instead of checking `chi_j(v_i)`, check that the **LCH14/HQC
butterfly** applied to a polynomial `f(x)` in the novelpoly basis produces
`(f(α+0), f(α+1), …, f(α+n-1))` and that pointwise multiplication of two such
forward FFTs followed by iFFT recovers the convolution.

## Recommended next steps

1. **Re-probe LCH14 with corrected multiplier.** Implement the HQC Algorithm 2
   butterfly verbatim using `s_{i-1}(a)` where `a` is a non-basis-element
   affine shift (e.g., `a = v_k` for `k > i`). Expected: 100% pass rate on
   the convolution-theorem probe. Sanity-check against Chen 2018 Algorithm 2
   line 9–10.
2. **Re-probe Gao-Mateer with the affine-subspace construction.** Use hamil
   Algorithm 3.5 instead of the Frobenius-DIF form. Expected: round-trip +
   convolution theorem both hold for `n ≤ 2^4` over GF(2^4).
3. **Port HQC Algorithm 2 to PAR3.** Target the canonical workload:
   `n = 4096` (slice+recovery), `m = 64`. Implementation:
   - BasisCvt (Algorithm 1 of Chen 2018): monomial → novelpoly via `s_i = x^2+x`-style recursion; XOR-only over GF(2^64) elements.
   - Butterfly: 1 PCLMULQDQ + 2 XORs per butterfly.
   - Inverse: same structure with `-= `flip (in char 2, `a → a+1`).
   - Datatype: `uint64_t` field elements with `x^64 + x^4 + x^3 + x + 1` irreducibles (matches Chen 2018 §4.3 and HQC §2.2).
4. **Benchmark vs PAR2.** Expected: PAR3-create reaches ~622 MB/s (matching
   PAR2), closing the 13.5× gap, on the canonical 1 GiB / 10K-slice / 1K-recovery.
5. **Update README and PHASE_2B_RESEARCH.md** to reflect that the open
   question is resolved, with citations to LCH14, Chen 2018, hamil, HQC,
   Polito.

## What this changes for parparpar issue #23

- Stage 0 (LCH14 verification) **can be reopened**, with the corrected probe
  expected to pass.
- Stage 1 (Gao-Mateer verification) **can be reopened**, with the affine
  variant expected to pass.
- The "no clean recursive FFT exists over GF(2^k)" verdict in
  PHASE_2B_RESEARCH.md "Root cause" section is reversed.
- Path A of the conclusion ("find or invent the Vandermonde sparse
  factorization") collapses to **Path A' = "implement HQC Algorithm 2
  verbatim, in C, with PCLMULQDQ"**, estimated 2–5 days, not 2–4 weeks.

## What I'm NOT claiming

- I have not re-run the probes myself to confirm 100% pass rate at `n=4`.
- I have not benchmarked HQC Algorithm 2 against the PAR3 kernel.
- Coxon 2021 (HAL 01845238) was blocked by Anubis; the LambdaClass blog was
  also blocked. The synthesis above relies on Chen 2018, hamil 2016, HQC 2026,
  and Polito 2024 as primary sources; LCH14 is cited via secondary references
  in those works.
- The HQC paper §2.3 quote uses field `F_{2^64}` with `m = 64` — exactly the
  PAR3 field. No field resizing required.

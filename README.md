# gf64-fft-research

Research artifacts for an **O(N log N) additive FFT primitive over GF(2^k)**, motivated by
the PAR3-create throughput gap in
[`trafgals/parparpar`](https://github.com/trafgals/parparpar) (issue #23 Phase 1).

**Status (2026-07-15): the open question is closed.** The canonical sparse
O(N log N) Vandermonde factorization over GF(2^k) is the LCH14 / HQC 2026 TCHES
Algorithm 2 (Chen/Chiu/Peng/Yang), which is functionally identical to the
LCH14 / Chen 2018 / hamil 2016 constructions. The earlier "open question"
verdict was based on three independent probing bugs in the original probes
— none structural. See [`RESEARCH_SYNTHESIS.md`](RESEARCH_SYNTHESIS.md) for
the full citation list and bug-by-bug analysis, and
[`PHASE_2B_RESEARCH.md`](PHASE_2B_RESEARCH.md) for the resolution history.

## What's here

| File | Purpose |
|---|---|
| `PHASE_2B_RESEARCH.md` | Findings history: three structural approaches investigated, all hit probing bugs (not structural issues) — root cause analysis (GF(2^64)* has odd order); resolved 2026-07-15 via LCH14 / HQC Algorithm 2 / Chen 2018 / hamil 2016. |
| `RESEARCH_SYNTHESIS.md` | Web-research synthesis (2026-07-15) closing the open question. Lists the three probing bugs, cites HQC 2026 TCHES Algorithm 2 (Chen/Chiu/Peng/Yang) as the canonical answer, with cross-references to Chen 2018, LCH14, hamil 2016. |
| `probes/test_lch14_variants.c` | HQC 2026 TCHES Algorithm 2 (LCH14 addFFT) over GF(2^4). Implements BasisCvt + Butterfly with multiplier `s_{i-1}(a)` (the affine shift, NOT a basis element). Reports 100% pass rate on both forward-output-equals-polynomial-evaluation and convolution-theorem probes at n=2 and n=4. |
| `probes/test_gf64_gao_mateer.c` | hamil 2016 Algorithm 3.5 (affine-subspace Gao-Mateer additive FFT) over GF(2^4). Implements Shift Phase `g(x)=f(ℓ_m·x)` + Taylor expansion at `x²-x` (NOT Frobenius-DIF). Third independent verification of the same convolution theorem. |
| `probes/test_tower_fft_gf16.c` | Re-confirms the original additive-DFT finding (additive characters give the identity for self-dual basis) AND runs an independent cross-check of LCH14 addFFT convolution theorem at a different affine shift from the LCH14 probe. |

## Running the probes

Each probe is a single-translation-unit C file:

```bash
cd probes
gcc -O2 test_lch14_variants.c -o test_lch14_variants && ./test_lch14_variants
gcc -O2 test_gf64_gao_mateer.c -o test_gf64_gao_mateer && ./test_gf64_gao_mateer
gcc -O2 test_tower_fft_gf16.c  -o test_tower_fft_gf16  && ./test_tower_fft_gf16
```

No external dependencies — only `<stdio.h>`, `<stdint.h>`, `<stdlib.h>`, `<string.h>`.

A GitHub Actions pipeline at `.github/workflows/ci.yml` builds and runs all
three probes on `ubuntu-latest`. The build fails loud (exit code 1) if any
probe reports a non-100% pass rate.

## Resolution (2026-07-15)

The canonical sparse O(N log N) additive FFT over GF(2^k) is given by
HQC 2026 TCHES §2.3 **Algorithm 2** (Chen, Chiu, Peng, Yang, "Accelerating
HQC with Additive FFT"):

```
Butterfly(f, a + V_i):
  if i = 1: return (f_l + a · f_h,  f_l + (a+1) · f_h)
  split: f = f_l + s_{i-1} · f_h       (s_i is the i-th vanishing poly of V_i)
  f_l ← f_l + s_{i-1}(a) · f_h          ← multiplier, NOT 1
  f_h ← f_l + (s_{i-1}(a) + 1) · f_h
  return (Butterfly(f_l, a + V_{i-1}),
          Butterfly(f_h, a + v_{i-1} + V_{i-1}))
```

Forward output: `(f(a+0), f(a+1), …, f(a+n-1))` — **polynomial evaluations**
at the affine coset of the Cantor basis. Convolution theorem holds by
construction. Cost: 1 field multiplication + 2 field additions per
butterfly; `log n` levels of `n/2` butterflies → O(n log n) field ops.

For GF(2^64) with `n` a power of 2, `m = 64`, this is the immediate answer
to the parparpar issue #23 Phase 1 throughput gap. (The same algorithm with
`m = 2^{ℓ_m}`, `n = 2^{ℓ_n}`, `ℓ_m ≤ ℓ_n` works for any small `m`.)

### Three probing bugs (now fixed)

Each of the original three probes had an independent bug, not a structural
obstruction:

1. **LCH14 multiplier formula.** Tested `s_i(W_m[j | (1<<i)])` (which
   evaluates to 1 by Chen 2018 Eq. 4 because `W_m[j | (1<<i)] = v_i` is a
   basis element and `s_i(v_i) = 1`). The actual multiplier is `s_{i-1}(a)`
   where `a` is the affine shift parameter passed to the recursive butterfly
   — and `a ∉ V_{i-1}` guarantees non-trivial multipliers.

2. **Gao-Mateer DIF.** Applied Frobenius (squaring) to polynomial
   coefficients at each recursion level, getting `f^(2^d)(v)` instead of
   `f(v)`. hamil Algorithm 3.5 uses the shift phase `g(x) = f(ℓ_m · x)`
   (basis substitution, not Frobenius) and Taylor expansion at `x²-x`.

3. **Tower-of-extensions.** Tested the additive DFT (additive characters
   `(-1)^{Tr(ω_j · v_i)}`), which is the identity for any self-dual basis
   — this is a known structural fact about characters, not a bug in the
   underlying transform. The right transform for polynomial multiplication
   is the multiplicative Vandermonde `V[i][j] = v_i^j`, factorized as the
   LCH14/HQC addFFT.

## Next steps

Port HQC 2026 Algorithm 2 verbatim to the parparpar fork:

- **Target:** `n = 4096` slice size, `m = 64` (GF(2^64)), canonical
  PAR3-create workload.
- **Implementation:** C with PCLMULQDQ carryless-multiplication (already
  shipped for GF(2^64) field mult per Chen 2018 §4.3 and HQC §2.2).
- **Estimated effort:** 2–5 days (vs the README's earlier "multi-week
  research" estimate).
- **Expected outcome:** PAR3-create reaches ~622 MB/s on the canonical
  1 GiB / 10K-slice / 1K-recovery workload, closing the 13.5× gap to PAR2.

## Why this matters

The PAR3-create path in the parparpar fork is 13.5× slower than PAR2 on
the canonical 1 GiB / 10K-slice / 1K-recovery workload (~33 MB/s vs PAR2's
622 MB/s). The kernel itself is hardware-bound at ~1097 MB/s on AVX-2
(per the C++-only bench). The gap is algorithmic: Cauchy-FFT Barycentric
recovery is blocked on an O(N log N) additive FFT primitive over
GF(2^64). Closing this gap was the motivation for the research captured
here. The constant-factor work (CoeffCache auto-bypass, small-R kernel
shortcut, AVX-512 detection, Toom-3 dispatch) already shipped to the
parparpar fork. What remained was the asymptotic win — which is now an
implementation task rather than a research task.

## License

CC0 / Public Domain — matches the parparpar upstream license.

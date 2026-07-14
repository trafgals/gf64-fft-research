# gf64-fft-research

Research artifacts for an **O(N log N) additive FFT primitive over GF(2^k)**, motivated by
the PAR3-create throughput gap in
[`trafgals/parparpar`](https://github.com/trafgals/parparpar) (issue #23 Phase 1).

The full investigation record is in [`PHASE_2B_RESEARCH.md`](PHASE_2B_RESEARCH.md). This repo
contains the standalone GF(2^4) probe programs that produced the findings; they do not depend
on the parparpar kernel and build/run with a single `gcc` command each.

## What's here

| File | Purpose |
|---|---|
| `PHASE_2B_RESEARCH.md` | Findings doc: three structural approaches (LCH14, Gao-Mateer DIF, tower-of-extensions) all hit fundamental issues in GF(2^4); root cause analysis (GF(2^64)* has odd order); Vandermonde sparse factorization identified as the open problem. |
| `probes/test_lch14_variants.c` | Tests three LCH14 multiplier variants (A/B/C) against a brute-force GF(2^4) convolution-theorem probe. All three give identical 12.10% pass rate (6975/57600) — the multiplier formula `s_i(W_m[j \| (1<<i)])` collapses to constant 1 in GF(2^4). |
| `probes/test_gf64_gao_mateer.c` | Decimation-in-frequency Gao-Mateer evaluation-basis FFT over GF(2^4). Round-trip works but the forward output is `f^(2^depth)(v)`, not polynomial evaluations — the "1/n scaling issue" prevents conv-theorem from holding in char 2. |
| `probes/test_tower_fft_gf16.c` | Enumerates evaluation sets V of size 4 in GF(2^4) and checks mixed trace + invertible pairing + self-dual basis. Finds 336 valid V's — but for every self-dual basis, the additive DFT matrix `F[i][j] = Tr(omega_j · V[i])` is the IDENTITY. The additive DFT IS the identity transform; it doesn't transform the input. |

## Running a probe

Each probe is a single-translation-unit C file:

```bash
cd probes
gcc -O2 test_lch14_variants.c -o test_lch14_variants && ./test_lch14_variants
gcc -O2 test_gf64_gao_mateer.c -o test_gf64_gao_mateer && ./test_gf64_gao_mateer
gcc -O2 test_tower_fft_gf16.c  -o test_tower_fft_gf16  && ./test_tower_fft_gf16
```

No external dependencies — only `<stdio.h>`, `<stdint.h>`, `<stdlib.h>`, `<string.h>`.

## The open question

The Vandermonde matrix `V[i][j] = v_i^j` (monomial basis evaluated at distinct `v_i`) IS the
correct transform for polynomial multiplication over GF(2^k):

- It IS invertible for distinct `v_i`.
- It DOES satisfy the convolution theorem: pointwise multiplication in the transformed basis
  = polynomial convolution in the monomial basis.

What's missing is its **sparse O(N log N) factorization** — the twiddle/multiplier structure
that lets a Cooley-Tukey-style recursive butterfly compute it. The literature does not
provide a closed form for this factorization in GF(2^k).

If you have a pointer to that factorization, or to a working implementation of an
`O(N log N)` additive FFT satisfying the convolution theorem over `GF(2^k)` for arbitrary k,
**this is exactly the missing piece** — please open an issue.

## Why this matters

The PAR3-create path in the parparpar fork is 13.5× slower than PAR2 on the canonical
1 GiB / 10K-slice / 1K-recovery workload (~33 MB/s vs PAR2's 622 MB/s). The kernel itself
is hardware-bound at ~1097 MB/s on AVX-2 (per the C++-only bench). The gap is algorithmic:
Cauchy-FFT Barycentric recovery is blocked on an O(N log N) additive FFT primitive over
GF(2^64).

Closing this gap is what motivated the research captured here. The constant-factor work
(CoeffCache auto-bypass, small-R kernel shortcut, AVX-512 detection, Toom-3 dispatch)
already shipped to the parparpar fork. What remains is the asymptotic win.

## License

CC0 / Public Domain — matches the parparpar upstream license.

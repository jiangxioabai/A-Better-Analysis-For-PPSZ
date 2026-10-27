# PPSZ numerical certificate

This package provides fixed-input, exact-rational checks for the numerical
estimates in *A Sharper Analysis of PPSZ for 3-SAT*. Its current entry point
is `checks/verify_repair.py`. It contains code and certificate data; the
manuscript and LaTeX sources are distributed separately.

## Reproduce

Use Python 3.10 or later. No third-party packages are required.

```sh
python checks/verify_repair.py
```

The expected first line is:

```text
PASS: 28 exact rational checks
```

The script prints its results and regenerates `checks/verification.json`.
The supplied `checks/verification_output.txt` is the genuine saved transcript
of the earlier successful run. The script, JSON and transcript are copied
unchanged from the current certificate; packaging did not rerun the checker.

## Files

| File | Purpose |
| --- | --- |
| `checks/verify_repair.py` | Fixed-input exact-rational verifier |
| `checks/verification.json` | Exact inputs, outward enclosures and passed-check names |
| `checks/verification_output.txt` | Saved actual execution output |
| `SHA256SUMS.txt` | Hashes of the four distributed payload files |

## Parameters and notation

The fixed inputs include `epsilon_R = 1/10`,
`epsilon_I = 0.07307238160252154687451293138`, `rho = 17/19`, and
`conditional_factor = 125/109`. The theorem bonuses remain
`gamma_target = 0.0000687793` and `eta_target = 0.000000364`, with safe
running-time bases `1.306969598` and `1.307031578`.

The JSON field `gamma_repair` contains the outward enclosure for the quantity
called `Gamma_avg` in the manuscript. Field names and the certificate
identifier `2026-10-27-structural-repair-integrated-v1` are retained for
compatibility. Generation metadata, including `publication_status`, are not
a live report of this repository's publication state.

## Scope

The checks cover finite arithmetic side conditions, coefficient signs,
entropy and cleanup constants, dual margins, the equality profile, fixed
lifting-root brackets and the stated base comparisons. The script performs
no parameter search or optimization. It is not a proof checker for the
combinatorial and probabilistic arguments, and no Lean formalization is
claimed.

The files under `checks/` form one matching certificate set. Earlier
certificates, if retained in a repository, should be identified as historical
and kept distinct from this entry point. The SHA-256 manifest describes the
bytes in this distribution; platform line-ending conversion changes text
hashes without changing the rational values.

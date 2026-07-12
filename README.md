# PPSZ manuscript revision bundle

## Files

- `ppsz_recombination_revised.tex` - complete arXiv-ready LaTeX source.
- `ppsz_recombination_revised.pdf` - compiled manuscript.
- `ppsz_certificate.json` - fixed exact-rational parameters, root brackets, and rounded targets.
- `verify_ppsz_constants.py` - independent exact-rational interval checker; performs no search or root finding.
- `verification_output.txt` - successful verifier transcript.
- `REVISION_NOTES.md` - substantive mathematical changes made in this revision.

## Verify the arithmetic

From this directory, run:

```sh
python3 verify_ppsz_constants.py ppsz_certificate.json
```

The final line must be:

```text
ALL CHECKS PASSED
```

The verifier uses `fractions.Fraction` for every proof-relevant endpoint. Logarithms are enclosed by a truncated atanh series with a rational remainder bound; exponentials are enclosed by a Taylor series with a rational geometric tail bound. Decimal conversion is used only for the human-readable transcript.

## Compile the manuscript

A standard pdfLaTeX installation is sufficient:

```sh
pdflatex -interaction=nonstopmode -halt-on-error ppsz_recombination_revised.tex
pdflatex -interaction=nonstopmode -halt-on-error ppsz_recombination_revised.tex
```

The source has no external bibliography, figures, or package-specific data files.

## Certified principal constants

- New Unique-3-SAT theorem gain: `0.0000687793`.
- Certified limiting recombination value: `0.0000687793804588365655...`.
- Scheder old unique gain: `1/15218 = 0.0000657116572479957944...`.
- Old lifted limiting gain: `0.0000003465837065168457...`.
- New lifted limiting gain: `0.0000003640269421506335...`.
- Relative increase in lifted gain: `5.032907002...%`.
- Safe general-3-SAT theorem gain: `0.000000364`.
- Safe outward-rounded general-3-SAT base: `1.307031578`.

## SHA-256 checksums

```text
17e74a0e297a3e8ec3d521d3272d3c85d2a44f1d60546d705ebc86811eb47982  ppsz_recombination_revised.tex
4ac9852b5e2d87cb18cbf39bfc346342e9c8310f1094b05f7e7a10814a2638ca  ppsz_recombination_revised.pdf
d683abf6fad7ed6b9983c782ae4308a66511cbfb938106b7dac1ee9ebb2aca5c  ppsz_certificate.json
1d29144b25b7a72de963787b587b804c443c3fef7ff3bfc33782a07ffd956be5  verify_ppsz_constants.py
1b5a8779d9d65bb785f895ec6c46fd4a94b0a987f3dca19ef6e88209c2ff0879  verification_output.txt
567c1f51345e3a436ad852eed4b80d83fe87b315d68e5702627239157bd6b86e  REVISION_NOTES.md
```

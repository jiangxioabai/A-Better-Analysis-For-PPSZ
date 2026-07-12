#!/usr/bin/env python3
"""Exact-rational interval verifier for the PPSZ recombination and lifting constants.

The script performs no parameter search and no root finding.  It reads fixed decimal
parameters and fixed rational brackets from ppsz_certificate.json, encloses logarithms
and exponentials by rational series with rigorous remainders, and checks every sign,
margin, and outward-rounded base used in the manuscript.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path


@dataclass(frozen=True)
class I:
    lo: Fraction
    hi: Fraction

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError("invalid interval")

    @staticmethod
    def point(x: Fraction | int) -> "I":
        q = x if isinstance(x, Fraction) else Fraction(x)
        return I(q, q)

    def __add__(self, other: "I | Fraction | int") -> "I":
        b = as_interval(other)
        return I(self.lo + b.lo, self.hi + b.hi)

    __radd__ = __add__

    def __neg__(self) -> "I":
        return I(-self.hi, -self.lo)

    def __sub__(self, other: "I | Fraction | int") -> "I":
        return self + (-as_interval(other))

    def __rsub__(self, other: "I | Fraction | int") -> "I":
        return as_interval(other) - self

    def __mul__(self, other: "I | Fraction | int") -> "I":
        b = as_interval(other)
        vals = (
            self.lo * b.lo,
            self.lo * b.hi,
            self.hi * b.lo,
            self.hi * b.hi,
        )
        return I(min(vals), max(vals))

    __rmul__ = __mul__

    def reciprocal(self) -> "I":
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError("interval contains zero")
        return I(Fraction(1, 1) / self.hi, Fraction(1, 1) / self.lo)

    def __truediv__(self, other: "I | Fraction | int") -> "I":
        return self * as_interval(other).reciprocal()

    def __rtruediv__(self, other: "I | Fraction | int") -> "I":
        return as_interval(other) / self


def as_interval(x: I | Fraction | int) -> I:
    if isinstance(x, I):
        return x
    return I.point(x)


def parse_q(s: str) -> Fraction:
    if "/" in s:
        a, b = s.split("/", 1)
        return Fraction(int(a), int(b))
    return Fraction(s)


def ln_unit_interval(y: Fraction, terms: int) -> I:
    """Enclose ln(y) for 1 <= y <= 2 by the atanh series."""
    if not (Fraction(1) <= y <= Fraction(2)):
        raise ValueError(f"ln_unit_interval argument out of range: {y}")
    if y == 1:
        return I.point(0)
    z = (y - 1) / (y + 1)
    z2 = z * z
    power = z
    partial = Fraction(0)
    for j in range(terms):
        partial += power / (2 * j + 1)
        power *= z2
    partial *= 2
    # power now equals z^(2N+1)
    tail = 2 * power / ((2 * terms + 1) * (1 - z2))
    return I(partial, partial + tail)


def scale_interval(a: I, k: int) -> I:
    return a * Fraction(k)


def ln_interval(x: Fraction, terms: int, ln2_cache: I | None = None) -> I:
    if x <= 0:
        raise ValueError("logarithm requires a positive argument")
    y = x
    k = 0
    while y < 1:
        y *= 2
        k -= 1
    while y > 2:
        y /= 2
        k += 1
    unit = ln_unit_interval(y, terms)
    ln2 = ln2_cache if ln2_cache is not None else ln_unit_interval(Fraction(2), terms)
    return unit + scale_interval(ln2, k)


def exp_point_interval(x: Fraction, terms: int) -> I:
    """Enclose exp(x) for 0 <= x < 1 by a Taylor series."""
    if not (0 <= x < 1):
        raise ValueError(f"exp input outside implemented range: {x}")
    term = Fraction(1)
    partial = Fraction(1)
    for j in range(1, terms + 1):
        term *= x
        term /= j
        partial += term
    next_term = term * x / (terms + 1)
    ratio = x / (terms + 2)
    tail = next_term / (1 - ratio)
    return I(partial, partial + tail)


def exp_interval(x: I, terms: int) -> I:
    lo = exp_point_interval(x.lo, terms).lo
    hi = exp_point_interval(x.hi, terms).hi
    return I(lo, hi)


def pow2_interval(x: I, ln2: I, exp_terms: int) -> I:
    return exp_interval(x * ln2, exp_terms)


def fkl_interval(t: Fraction, log_terms: int, ln2: I) -> I:
    return (1 - t) * ln_interval(1 - t, log_terms, ln2) + t


def h2_interval(d: Fraction, log_terms: int, ln2: I) -> I:
    if not (0 < d < 1):
        raise ValueError("binary entropy argument must lie in (0,1)")
    nat = -(d * ln_interval(d, log_terms, ln2) + (1 - d) * ln_interval(1 - d, log_terms, ln2))
    return nat / ln2


def format_fraction(q: Fraction, significant_digits: int = 22) -> str:
    # Reporting only.  Decimal conversion is performed at high precision and never
    # participates in a proof check.
    with localcontext() as ctx:
        ctx.prec = significant_digits + 12
        value = Decimal(q.numerator) / Decimal(q.denominator)
        return format(value, f".{significant_digits}g")


def sci_interval(a: I, places: int = 22) -> str:
    return f"[{format_fraction(a.lo, places)}, {format_fraction(a.hi, places)}]"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def require_positive(a: I, message: str) -> None:
    require(a.lo > 0, f"{message}: lower endpoint is not positive")


def require_negative(a: I, message: str) -> None:
    require(a.hi < 0, f"{message}: upper endpoint is not negative")


def main() -> int:
    root = Path(__file__).resolve().parent
    if len(sys.argv) > 2:
        raise ValueError("usage: verify_ppsz_constants.py [certificate.json]")
    cert_path = Path(sys.argv[1]).resolve() if len(sys.argv) == 2 else root / "ppsz_certificate.json"
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    log_terms = int(cert["series"]["log_terms"])
    exp_terms = int(cert["series"]["exp_terms"])

    ln2 = ln_unit_interval(Fraction(2), log_terms)
    p0 = 2 * ln2 - 1
    pstar = 1 - 1 / (2 * ln2)
    q0 = p0 - pstar
    require_positive(q0, "q0 = p0-pstar")

    # Main Unique-3-SAT recombination.
    u = cert["unique"]
    eps_r = parse_q(u["epsilon_R"])
    eps_i = parse_q(u["epsilon_I"])
    gamma_reported = parse_q(u["gamma_reported"])
    require(Fraction(0) < eps_r < Fraction(13, 100), "epsilon_R range")
    require(Fraction(0) < eps_i < Fraction(1, 5), "epsilon_I range")

    c_l = parse_q("0.001687") * eps_r - parse_q("0.006404") * eps_r * eps_r
    A = Fraction(17, 18) * c_l
    threshold = 2 * A / parse_q("0.9")
    p_reg = parse_q("1.1") * eps_r * threshold
    c_t = I.point(parse_q("0.009307") - parse_q("0.055") * eps_r) - parse_q("0.1503") * fkl_interval(eps_r, log_terms, ln2)
    S = c_t - 5 * A

    f_i = fkl_interval(eps_i, log_terms, ln2)
    f_5i = fkl_interval(5 * eps_i, log_terms, ln2)
    b1 = I.point(parse_q("0.030966") * eps_i - parse_q("0.0028") * eps_i * eps_i) - parse_q("0.4027") * f_i
    b0 = I.point(parse_q("0.06259") * eps_i) - parse_q("0.344") * f_i
    bT = I.point(parse_q("0.009307") - parse_q("0.2405") * eps_i - parse_q("0.03125") * eps_i * eps_i) - parse_q("0.06183") * f_5i

    require(A > p_reg > 0, "A>P_reg>0")
    require_positive(S, "S>0")
    require_positive(b1, "b1>0")
    require_positive(b0, "b0>0")
    require_negative(bT, "bT<0")
    dual_i0 = b0 - 2 * b1
    dual_tau_cross = A * bT + b1 * S
    require_positive(dual_i0, "b0-2b1>0")
    require_positive(dual_tau_cross, "A*bT+b1*S>0")

    lam = b1 / A
    dual_tau = bT + lam * S
    require_positive(dual_tau, "bT+lambda*S>0")
    gamma_star = b1 * (A - p_reg) / (A + b1)
    require(gamma_star.lo > gamma_reported, "gamma_star exceeds reported gamma")
    tight_i1 = (A - p_reg) / (A + b1)
    require(tight_i1.hi < 1, "tight witness lies in the simplex")

    # Enclosures printed in the manuscript.  These are deliberately wider than
    # the internal intervals, and their directions are checked explicitly.
    rep = cert["reported_intervals"]
    def contained(value: I, lower_key: str, upper_key: str, label: str) -> None:
        require(parse_q(rep[lower_key]) <= value.lo, f"{label} displayed lower endpoint")
        require(value.hi <= parse_q(rep[upper_key]), f"{label} displayed upper endpoint")

    contained(I.point(A), "A_lower", "A_upper", "A")
    contained(I.point(p_reg), "P_reg_lower", "P_reg_upper", "P_reg")
    contained(S, "S_lower", "S_upper", "S")
    contained(b1, "b1_lower", "b1_upper", "b1")
    contained(b0, "b0_lower", "b0_upper", "b0")
    contained(bT, "bT_lower", "bT_upper", "bT")
    contained(lam, "lambda_lower", "lambda_upper", "lambda")
    contained(gamma_star, "gamma_star_lower", "gamma_star_upper", "gamma_star")
    contained(tight_i1, "tight_i1_lower", "tight_i1_upper", "tight_i1")

    unique_base = pow2_interval(p0 - gamma_reported, ln2, exp_terms)
    require(unique_base.hi < parse_q(u["unique_base_upper"]), "outward unique-case base")

    # Quantitative lifting function.
    l = cert["lifting"]
    gamma_old = parse_q(l["gamma_old"])
    gamma_new = parse_q(l["gamma_new"])
    require(gamma_new > gamma_old, "new unique gain strictly exceeds Scheder's gain")
    old_unique_base = pow2_interval(p0 - gamma_old, ln2, exp_terms)
    require(
        old_unique_base.hi < parse_q(l["old_unique_base_upper"]),
        "Scheder unique-case base outward rounding",
    )
    unique_relative = gamma_new / gamma_old - 1
    require(
        unique_relative > parse_q(l["unique_relative_improvement_lower"]),
        "relative unique-gain improvement",
    )

    def root_function(d: Fraction, gamma: Fraction) -> I:
        # h_2(d) + (1-p_*+gamma)d - gamma
        return h2_interval(d, log_terms, ln2) + (1 - pstar + gamma) * d - gamma

    def certify_lift(label: str, gamma: Fraction, dlo_s: str, dhi_s: str) -> tuple[I, I, I]:
        dlo = parse_q(dlo_s)
        dhi = parse_q(dhi_s)
        require(Fraction(0) < dlo < dhi < Fraction(1, 2), f"{label} root bracket")
        require_negative(root_function(dlo, gamma), f"{label} root lower sign")
        require_positive(root_function(dhi, gamma), f"{label} root upper sign")
        delta = I(dlo, dhi)
        eta = q0 * delta
        base = pow2_interval(p0 - eta, ln2, exp_terms)
        return delta, eta, base

    delta_old, eta_old, base_old = certify_lift(
        "old", gamma_old, l["delta_old_lower"], l["delta_old_upper"]
    )
    delta_new, eta_new, base_new = certify_lift(
        "new", gamma_new, l["delta_new_lower"], l["delta_new_upper"]
    )
    require(eta_new.lo > eta_old.hi, "new lifted gain strictly exceeds old lifted gain")

    relative = eta_new / eta_old - 1
    require(relative.lo > parse_q(l["relative_improvement_lower"]), "relative lifted-gain improvement")
    require(base_old.hi < parse_q(l["old_base_upper"]), "old lifted base outward rounding")
    require(base_new.hi < parse_q(l["new_limiting_base_upper"]), "new limiting base outward rounding")

    # Safe theorem constant and the two limiting branch margins at a fixed rational delta.
    delta_safe = parse_q(l["delta_safe"])
    eta_safe = parse_q(l["eta_safe"])
    high_branch = q0 * delta_safe
    low_branch = gamma_new * (1 - delta_safe) - (1 - p0) * delta_safe - h2_interval(delta_safe, log_terms, ln2)
    high_margin = high_branch - eta_safe
    low_margin = low_branch - eta_safe
    require_positive(high_margin, "safe high-I branch margin")
    require_positive(low_margin, "safe unique-residual branch margin")
    claimed_margin = parse_q(rep["safe_branch_margin_lower"])
    require(high_margin.lo > claimed_margin, "reported high-branch safety margin")
    require(low_margin.lo > claimed_margin, "reported low-branch safety margin")
    safe_base = pow2_interval(p0 - eta_safe, ln2, exp_terms)
    require(safe_base.hi < parse_q(l["safe_general_base_upper"]), "safe theorem base outward rounding")

    # Fixed decimal enclosures quoted in the paper.
    require(Fraction("0.0000003465837065") <= eta_old.lo, "old eta lower displayed enclosure")
    require(eta_old.hi <= Fraction("0.0000003465837066"), "old eta upper displayed enclosure")
    require(Fraction("0.0000003640269421") <= eta_new.lo, "new eta lower displayed enclosure")
    require(eta_new.hi <= Fraction("0.0000003640269422"), "new eta upper displayed enclosure")
    require(relative.lo > Fraction("0.0503290"), "5.03290 percent lower enclosure")

    print(f"certificate version: {cert['version']}")
    print(f"ln 2                       {sci_interval(ln2, 20)}")
    print(f"p0=2 ln 2-1               {sci_interval(p0, 20)}")
    print(f"p*=1-1/(2 ln 2)           {sci_interval(pstar, 20)}")
    print(f"q0=p0-p*                  {sci_interval(q0, 20)}")
    print()
    print("Unique-3-SAT recombination")
    print(f"A                          {sci_interval(I.point(A), 20)}")
    print(f"P_reg                      {sci_interval(I.point(p_reg), 20)}")
    print(f"S                          {sci_interval(S, 20)}")
    print(f"b1                         {sci_interval(b1, 20)}")
    print(f"b0                         {sci_interval(b0, 20)}")
    print(f"bT                         {sci_interval(bT, 20)}")
    print(f"b0-2b1                     {sci_interval(dual_i0, 20)}")
    print(f"bT+lambda S                {sci_interval(dual_tau, 20)}")
    print(f"gamma_*                    {sci_interval(gamma_star, 20)}")
    print(f"gamma_*-gamma_reported     {sci_interval(gamma_star-gamma_reported, 20)}")
    print(f"unique base                {sci_interval(unique_base, 20)}")
    print(f"Scheder unique gain        {format_fraction(gamma_old, 20)}")
    print(f"Scheder unique base        {sci_interval(old_unique_base, 20)}")
    print(f"unique-gain improvement    {format_fraction(unique_relative, 20)}")
    print()
    print("Scheder-Steinberger lifting")
    print(f"delta_old                  {sci_interval(delta_old, 20)}")
    print(f"eta_old                    {sci_interval(eta_old, 20)}")
    print(f"old general base           {sci_interval(base_old, 20)}")
    print(f"delta_new                  {sci_interval(delta_new, 20)}")
    print(f"eta_new                    {sci_interval(eta_new, 20)}")
    print(f"new limiting general base {sci_interval(base_new, 20)}")
    print(f"eta_new/eta_old-1          {sci_interval(relative, 20)}")
    print()
    print("Safe general theorem")
    print(f"delta_safe                 {delta_safe}")
    print(f"high-branch margin         {sci_interval(high_margin, 20)}")
    print(f"low-branch margin          {sci_interval(low_margin, 20)}")
    print(f"safe general base          {sci_interval(safe_base, 20)}")
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # concise failure for CI and artifact checking
        print(f"VERIFICATION FAILED: {exc}", file=sys.stderr)
        raise

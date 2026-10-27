#!/usr/bin/env python3
"""Exact-rational verification of the PPSZ structural-repair constants.

This checks finite numerical claims, not the combinatorial/probabilistic proof.
No floating point is used in any decision. Only the Python standard library is
required. Run: python verify_repair.py
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as Q
from math import comb, factorial
from pathlib import Path
import json

TERMS = 90

@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q
    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError('Reversed interval')
    @staticmethod
    def lift(x):
        if isinstance(x, Interval):
            return x
        return Interval(Q(x), Q(x))
    def __add__(self, other):
        other = self.lift(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)
    __radd__ = __add__
    def __neg__(self):
        return Interval(-self.hi, -self.lo)
    def __sub__(self, other):
        return self + (-self.lift(other))
    def __rsub__(self, other):
        return self.lift(other) - self
    def __mul__(self, other):
        other = self.lift(other)
        ends = [x*y for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        return Interval(min(ends), max(ends))
    __rmul__ = __mul__
    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError('Interval includes zero')
        return Interval(1/self.hi, 1/self.lo)
    def __truediv__(self, other):
        return self * self.lift(other).reciprocal()
    def __rtruediv__(self, other):
        return self.lift(other) / self

def unit_log(x: Q) -> Interval:
    """log(x) for rational x in [1,2], using an atanh series and tail."""
    if not 1 <= x <= 2:
        raise ValueError('unit_log requires 1 <= x <= 2')
    z = (x-1)/(x+1)
    s = 2*sum((z**(2*j+1)/Q(2*j+1) for j in range(TERMS)), Q(0))
    tail = 2*z**(2*TERMS+1)/(Q(2*TERMS+1)*(1-z*z))
    return Interval(s, s+tail)

LN2 = unit_log(Q(2))

def logq(x: Q) -> Interval:
    if x <= 0:
        raise ValueError('log argument must be positive')
    y, k = x, 0
    while y < 1:
        y *= 2
        k -= 1
    while y > 2:
        y /= 2
        k += 1
    return unit_log(y) + k*LN2

def fkl(x: Q) -> Interval:
    if not 0 <= x < 1:
        raise ValueError('fkl requires 0 <= x < 1')
    return (1-x)*logq(1-x) + x

def expq(x: Q) -> Interval:
    if x < 0 or x >= TERMS+1:
        raise ValueError('expq requires 0 <= x < TERMS+1')
    s = sum((x**j/Q(factorial(j)) for j in range(TERMS)), Q(0))
    first = x**TERMS/Q(factorial(TERMS))
    tail = first/(1-x/Q(TERMS+1))
    return Interval(s, s+tail)

def expi(x: Interval) -> Interval:
    return Interval(expq(x.lo).lo, expq(x.hi).hi)

def outward(x, digits=23):
    x = Interval.lift(x)
    scale = 10**digits
    a = x.lo.numerator*scale // x.lo.denominator
    b = -((-x.hi.numerator*scale) // x.hi.denominator)
    def fmt(n):
        sign = '-' if n < 0 else ''
        n = abs(n)
        return f'{sign}{n//scale}.{n%scale:0{digits}d}'
    return [fmt(a), fmt(b)]

checks = []
def require(condition, name):
    if not condition:
        raise AssertionError(name)
    checks.append(name)

# The series tail uses that (k/(k+2))**k decreases with k.
# For d >= 101, k=d+1 >= 102, so its value is below 1/7.
require(7*102**102 < 104**102, 'series tail ratio < 1/7 at k=102')
series_part = sum((Q(d+1,d+3)**(d+1)/Q((d+3)**2)
                   for d in range(1,101)), Q(0))
series_upper = series_part + Q(1,7*102)
require(series_upper < Q(109,2000), 'sum_d (d+1)^(d+1)/(d+3)^(d+3) < 0.0545')

eR = Q(1,10)
eI = Q('0.07307238160252154687451293138')
rho = Q(17,19)
K = Q(1)/(1-Q(2,25)*eR*16)
require(K == Q(125,109), 'conditional amplification = 125/109')
# phi in [-1/sqrt(5),1]; 17*eR < sqrt(5) ensures density positivity.
require((17*eR)**2 < 5, 'linear component densities are strictly positive')
# max_{r in [1/5,1/2]}(5r-1)(1-2r)=9/40, by completing a square.
loss_factor = Q(51,50)*K*(Q(9,20)+8*Q(109,2000))
C = Q(1037,1000)
require(loss_factor < C, 'refined cleanup loss constant < 1.037')

# Exact moments of phi, using t=sqrt(1-2r).
def moment(d):
    return sum((Q(comb(d,j))*5**j*(-3)**(d-j)/Q(2**d*(d+2*j+2))
                for j in range(d+1)), Q(0))
m1,m2,m3,m4 = [moment(d) for d in (1,2,3,4)]
require((m1,m2,m3,m4) == (Q(0),Q(3,32),Q(12,385),Q(9,224)), 'moments of phi')
quartic_common = m4*m4+6*m4*m2*m2+42*m2**4
quartic_triangle = m4*m4+6*m4*m2*m2+12*m3*m3*m2
quartic_fourcycle = m4*m4+6*m4*m2*m2+9*m2**4
require(quartic_triangle <= quartic_common, 'triangle fourth moment covered')
require(quartic_fourcycle <= quartic_common, '4-cycle fourth moment covered')
kl_numerator = m2*m2/2 + eR*eR*quartic_common/3
require(LN2.lo > Q(6931,10000), 'log(2) > 0.6931')
require(kl_numerator < Q('0.006374')*Q('0.6931'), 'KL per edge < 0.006374 epsilon_R^2')

cL = Q('0.001687')*eR - Q('0.006378')*eR*eR
A = rho*cL
Thr = 2*A/Q('0.9')
P = C*rho*eR*Thr
cT = Q('0.009307')-Q('0.055')*eR-Q('0.1503')*fkl(eR)
# Conservative global TwoCC-free convention on multigraphs: 4|TwoCC|,
# plus 2|TwoCC| from converting J_i into ID_i.
S = cT-6*A
b1 = Q('0.030966')*eI-Q('0.0028')*eI*eI-Q('0.4027')*fkl(eI)
b0 = Q('0.06259')*eI-Q('0.344')*fkl(eI)
bT = Q('0.009307')-Q('0.2405')*eI-Q('0.03125')*eI*eI-Q('0.06183')*fkl(5*eI)
lam = b1/A
gamma = b1*(A-P)/(A+b1)
gamma_target = Q('0.0000687793')
require(0 < eR <= Q('0.1'), 'regular parameter inside printed 0.1 range')
require(0 < eI <= Q(1,5), 'irregular density parameter inside admissible range')
require(0 < Thr < Q(1,1150), 'high-density threshold condition')
require(A>P>0, 'A > P > 0')
require(S.lo>0 and b1.lo>0 and b0.lo>0 and bT.hi<0, 'coefficient signs')
require((b0-2*b1).lo>0, 'dual i0 slack > 0')
require((bT+lam*S).lo>0, 'dual TwoCC slack > 0')
require(gamma.lo > gamma_target+Q('0.000000096'), 'repaired gain exceeds target by more than 9.6e-8')
require(gamma_target > Q(1,15218), 'original target still exceeds published bonus')

# Unchanged target-to-base / lifting comparisons, independently enclosed.
p0=2*LN2-1
pstar=1-1/(2*LN2)
q0=p0-pstar
eta=Q('0.000000364')
delta=Q('0.00000338183369')
h2=(-delta*logq(delta)-(1-delta)*logq(1-delta))/LN2
ubranch=gamma_target*(1-delta)-(1-p0)*delta-h2
highbranch=q0*delta
require((highbranch-eta).lo>Q('0.0000000000269'), 'lifting high-branch margin > 2.69e-11')
require((ubranch-eta).lo>Q('0.0000000000269'), 'lifting unique-residual margin > 2.69e-11')
base_u=expi((p0-gamma_target)*LN2)
base_g=expi((p0-eta)*LN2)
require(base_u.hi<Q('1.306969598'), 'unique target base < 1.306969598')
require(base_g.hi<Q('1.307031578'), 'general target base < 1.307031578')


# Integration-only additions: the unchanged lifting-root inputs and the new
# equality profile. These are fixed-input checks, not a search or optimization.
equality_i1=(A-P)/(A+b1)
require(equality_i1.lo>0 and equality_i1.hi<1, 'relaxed equality profile lies in (0,1)')

def entropy_at(x):
    return (-x*logq(x)-(1-x)*logq(1-x))/LN2

def lifting_root_expression(x,gain):
    return entropy_at(x)+(1-pstar+gain)*x-gain

root_inputs={
    'old':(Q(1,15218),Q('0.00000321978491531273261'),Q('0.00000321978491531273262')),
    'target':(gamma_target,Q('0.00000338183369577144614'),Q('0.00000338183369577144615')),
}
root_results={}
for name,(gain,lower,upper) in root_inputs.items():
    require(lifting_root_expression(lower,gain).hi<0, f'{name} lifting root lower endpoint has negative sign')
    require(lifting_root_expression(upper,gain).lo>0, f'{name} lifting root upper endpoint has positive sign')
    lifted=q0*Interval(lower,upper)
    root_results['eta_limit_'+name]=lifted
    root_results['limiting_base_'+name]=expi((p0-lifted)*LN2)

values = dict(epsilon_R=eR,epsilon_I=eI,rho=rho,conditional_factor=K,
              cleanup_bound=loss_factor,cleanup_used=C,series_upper=series_upper,
              KL_coefficient_upper=kl_numerator/Q('0.6931'),
              cL=cL,A=A,threshold=Thr,penalty=P,cT=cT,S=S,b0=b0,b1=b1,bT=bT,
              dual_weight=lam,dual_i0_slack=b0-2*b1,dual_TwoCC_slack=bT+lam*S,
              gamma_repair=gamma,gamma_target=gamma_target,gain_slack=gamma-gamma_target,
              unique_base=base_u,general_base=base_g,
              lifting_high_slack=highbranch-eta,lifting_residual_slack=ubranch-eta)
values.update(equality_i1=equality_i1,A_minus_P=A-P,q0=q0,
              eta_target=eta,separator=delta,**root_results)
result={
    'certificate_version':'2026-10-27-structural-repair-integrated-v1',
    'publication_status':'local integration bundle; not claimed published to GitHub',
    'exact_inputs':{
        'epsilon_R':'1/10',
        'epsilon_I':'0.07307238160252154687451293138',
        'rho':'17/19',
        'conditional_factor':'125/109',
        'cleanup_used':'1.037',
        'entropy_bound':'0.006374',
        'entropy_used':'0.006378',
        'gamma_target':'0.0000687793',
        'eta_target':'0.000000364',
        'separator':'0.00000338183369',
        'delta_old_lower':'0.00000321978491531273261',
        'delta_old_upper':'0.00000321978491531273262',
        'delta_target_lower':'0.00000338183369577144614',
        'delta_target_upper':'0.00000338183369577144615'
    },
    'status':'PASS: exact rational arithmetic checks only; not a proof checker',
    'series_terms':TERMS,
    'checks':checks,
    'outward_decimal_enclosures':{k:outward(v) for k,v in values.items()}
}
output = Path(__file__).with_name('verification.json')
output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(f'PASS: {len(checks)} exact rational checks')
for k in ['rho','cleanup_bound','KL_coefficient_upper','A','threshold','penalty','S',
          'dual_i0_slack','dual_TwoCC_slack','gamma_repair','gain_slack',
          'unique_base','general_base','dual_weight','equality_i1',
          'eta_limit_old','eta_limit_target','limiting_base_old','limiting_base_target']:
    lo,hi=outward(values[k])
    print(f'{k}: [{lo}, {hi}]')
print('Proof coverage: the manuscript supplies the combinatorial/probabilistic arguments; this transcript checks arithmetic only.')

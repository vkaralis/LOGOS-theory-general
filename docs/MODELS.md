# Mathematical models and comparison contracts

All examples are dimensionless. Bolus administrations are instantaneous dose
atoms; a continuous-rate bound does not bound an ideal atom.

## Finite Emax model

Between atoms, dC/dt = u/V - k C; C(t_i+) = C(t_i-) + D_i/V.
Benefit b(C) = Emax C/(EC50+C), J = integral_0^T b(C(t)) dt.
Set V=k=EC50=Emax=D=1, T=3 and C(0-)=0.

Classes: equal boluses at 0,1,2; optimized amounts at these times; amounts and
ordered times (first at zero); loading plus bounded continuous input and coast.
All have total administered amount one, C <= 1, u <= 1 and injury
z'=C²-0.5z from z(0)=0, z <= 0.5 over the whole future tail.
The shared mass-kernel bound gives C <= 1 and integral_0^infinity C² dt <= 1/2,
so z <= 1/2 for every nonnegative compared schedule, beyond any plotting horizon.

For exponential segments J increments by log((1+Cstart)/(1+Cend)).
SLSQP uses ftol=1e-13 and maxiter=1000. Timing uses 24 starts, seed 20261003.
Successful optimization and feasibility are checked independently.

The expanded design loads Vc, maintains u=kVc for L, then coasts for delta:

- delta = log(2+EC50/c)/k
- L = T-delta
- D = Vc[1+kT-log(2+EC50/c)].

For T <= log(2+EC50 V/D)/k the optimum is a single initial atom.
Concavity and a marginal support inequality certify the ideal solution. Timing
optimization remains a multistart numerical candidate. Expected J values are
0.629606362035524, 0.672945396127757, 0.677043273327298, 0.679552690590195.

## Periodic sensitivity and injury

Mean administered input is 0.5; k=V=1. A cycle of length tau administers 0.5 tau,
with C(t)=0.5 tau exp(-t)/(1-exp(-tau)). State equations:

- a'=r(C)(1-a)-d(C)a; benefit=a q(C).
- z'=v(C)-rho(C)z.
- q_conc=C/(1+C); q_sig=C³/(1+C³).
- d(C)=2q(C); v(C)=C⁸/(0.35⁸+C⁸).
- Independent: r=r0, rho=rho0.
- Dependent: r=r0(1-q)^4, rho=rho0(1-v).

Baseline r0=0.7, rho0=1, peak C <= 4, peak z <= 3.5.
Cases A–D are concave/independent, sigmoidal/independent,
concave/dependent and sigmoidal/dependent. All use identical injury generation.

Affine integrating factors determine each periodic initial state. Simpson
quadrature uses 8193 points for searches and 40001 at reference intervals;
selected results also undergo an independent augmented-ODE calculation.
Search: 121 intervals in [0.08,12], sign-changing cap boundaries refined by
brentq and objective refinement near the best detected feasible candidate.
This reconstruction can miss isolated feasible components or additional maxima.
It does not reproduce an unavailable historical optimization implementation.

The mean-benefit negative control uses case A:
mean(a q) <= r0 q(mean C)/(r0+2q(mean C)), attained by constant exposure.
The other cases alter the conditions of this bound. Constant exposure can be
infeasible, so an efficacy ratio against that exposure is not a gain against a
feasible comparator. Periodic closure excludes startup and shutdown transients.

The local grid uses r0 in {0.56,0.70,0.84}, rho0 in {0.8,1.0,1.2} and cap in
{2.8,3.5,4.2}; additional caps 1 and 20 are 18 stress tests. Retain all outcomes,
including no feasible cycle detected at cap 1.

## Cessation tail

C(t)=C0 exp(-kt); z'=alpha C²-rho z, with C0,z0,alpha,rho>=0 and k>0.
For rho != 2k, A=alpha C0²/(2k-rho),
z(t)=(z0+A)exp(-rho t)-A exp(-2kt).
For rho=2k, z(t)=(z0+alpha C0²t)exp(-2kt).
The peak compares time zero, every admissible positive stationary time and the
limiting endpoint. Zero repair yields an asymptotic supremum.

Example C0=alpha=k=1, rho=0.2, z0=0: injury at 0.05 is 0.047340, while the
future peak is 0.387132 at 1.279214. The cap 0.1 is violated after cessation.

## Additional verification models

`independent_benchmarks.json` also separates quadratic finite-allocation,
linear injury and conditional-wait analytical examples. Their parameters,
resource rules and objectives are defined in `benchmarks.py`; they must not be
combined with the fixed-amount Emax or periodic results. These are mathematical
implementation checks only. No clinical prediction metrics or outcome claims
are generated from the synthetic cohort.

## Relaxed Hill envelope

For b(C)=Emax C^h/(EC50^h+C^h), h>1, the tangent concentration is
Cdagger=EC50(h-1)^(1/h). On [0,Cmax], support=min(Cmax,Cdagger).
The envelope is the zero-to-support chord below support and the Hill curve
above support. `static.py` provides the envelope; six test cases check dominance,
concavity and supporting slopes. This is a static relaxation, not a physically
attainable PK trajectory or a regimen optimization.

## Recovery-gated rectangular inputs (separate model)

`recovery_models.py`, `recovery.py` and `recovery_finite.py` use C'=kappa(u-C),
E=C/(1+C), R'=r W(C)(1-R)-alpha E R,
Z'=beta E-pi Wz(C) Z. The outcome rate is eta E R-rho(1-E).
Here rho penalizes untreated burden; it is not the injury-repair rate.
W=1/[1+(C/chi)^n] and Wz=1/[1+(C/chiZ)^2]. Parameters are
kappa=4, alpha=2, pi=0.7, beta=0.2, eta=1, rho=0.08, chi=0.35,
chiZ=1, n=4 and r=1.

Stationary inputs use average-input cap 1, input cap 5 and injury cap 0.22.
The rectangular witness uses period 1.679903, duty 0.2, on-input 5, off-input 0.
A one-million-interval derivative grid bounds the entire admissible constant
family; recovery-clamped and recovery-disabled negative controls are retained.
A deterministic 32-point Sobol parameter box (±20%) and 63 candidate rectangles
per vector retain negative/unsafe outcomes as well as positive outcomes.

The finite comparison starts C=0,R=1,Z=0, permits administration to time 20 and
evaluates outcome through 25. Twelve on/off cycles use duty 0.2. The whole
constant-input family u in [0,1] is bounded with 4096 intervals and a derivative
allowance; independent DOP853/Radau integrations and declining-tail conditions
verify the witness. Common available amount cap is 20; realized amounts differ.
These are floating-point bounds, not interval-arithmetic certificates.

Imperfect information uses three equally likely recovery rates {0.5,1,2},
wait bounds [0,6], utility 0.8(1-exp(-r w))-0.15(exp(0.4w)-1), and hypothetical
log-rate error SD {0,0.25,0.5,1,2}. Gauss–Hermite quadrature at 80/160 nodes
checks convergence. Net values at costs 0.002 and 0.01 retain unfavorable cases.
No empirical sensor accuracy, clinical policy effect or prediction is estimated.

## Fixed-target minimum-dose frontier (separate model)

`frontier.py`: C'=u-C, Z'=C-0.2 Z, benefit 4C-C²/2 integrated through time 4.
Z starts at zero; C can have an initial nonnegative atom. Z <= 2.5 over the
entire future path. Flexible input has an initial atom and 80 piecewise-constant
bins, with rates <=20. The search enforces a 1e-6 injury slack, minimizes total
amount at targets {10,11,12,13}, and refines sign-changing peak brackets.
Independent DOP853 integration checks benefit and injury through time 20;
zero input and declining analytic tails exclude later maxima.

The constant comparator starts without an atom, has the same benefit horizon
and full-tail cap, and is optimized analytically over all constant rates.
Its maximum benefit is 11.270930. Targets 10 and 11 show about 20% dose saving;
12 and 13 are infeasible for that constant family, so percentages are undefined.
Concavity supports an LP lower bound within the declared finite input family.
The bounds do not prove an unrestricted continuous-input minimum.

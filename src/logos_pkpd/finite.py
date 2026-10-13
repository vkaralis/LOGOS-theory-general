"""Finite cumulative Emax benefit at a fixed administered amount."""
import numpy as np
from scipy.optimize import minimize, brentq

SEED = 20261003

def objective(doses, times, horizon=3.0):
    doses, times = np.asarray(doses), np.asarray(times)
    if len(doses) != len(times) or np.any(doses < 0) or np.any(np.diff(times) < 0) or np.any(times < 0) or np.any(times > horizon):
        raise ValueError('Doses must be nonnegative and times ordered inside the horizon')
    return _integral(doses, times, horizon)

def _integral(doses, times, horizon):
    """Smooth segment integral used at constrained solver trial points."""
    c = value = previous = 0.0
    for dose, time in zip(doses, times):
        value += np.log1p(c) - np.log1p(c * np.exp(-(time - previous)))
        c = c * np.exp(-(time - previous)) + dose
        previous = time
    return float(value + np.log1p(c) - np.log1p(c * np.exp(-(horizon - previous))))

def plateau(dose=1., horizon=3., k=1., volume=1., ec50=1.):
    if min(dose, horizon, k, volume, ec50) <= 0:
        raise ValueError('Physical parameters must be positive')
    threshold = np.log(2 + ec50 * volume / dose) / k
    if horizon <= threshold:
        return dict(concentration=dose / volume, maintenance=0., coast=horizon)
    c = brentq(lambda c: volume*c*(1+k*horizon-np.log(2+ec50/c))-dose,
               dose/(volume*(1+k*horizon)), dose/volume, xtol=1e-14)
    coast = np.log(2+ec50/c)/k
    return dict(concentration=float(c), maintenance=float(horizon-coast), coast=float(coast))

def run():
    options = dict(ftol=1e-13, maxiter=1000)
    fixed = minimize(lambda x: -objective(x, [0,1,2]), np.ones(3)/3,
                     method='SLSQP', bounds=[(0,1)]*3,
                     constraints={'type':'eq','fun':lambda x: x.sum()-1}, options=options)
    if not fixed.success:
        raise RuntimeError(fixed.message)
    rng = np.random.default_rng(SEED)
    starts = [np.r_[fixed.x,1.,2.]] + [np.r_[rng.dirichlet(np.ones(3)), np.sort(rng.uniform(0,3,2))] for _ in range(23)]
    candidates = []
    for start in starts:
        fit = minimize(lambda x: -_integral(x[:3], [0,*x[3:]], 3.),
                       start, method='SLSQP', bounds=[(0,1)]*3+[(0,3)]*2,
                       constraints=[{'type':'eq','fun':lambda x: x[:3].sum()-1},
                                    {'type':'ineq','fun':lambda x: x[4]-x[3]}], options=options)
        if fit.success and abs(fit.x[:3].sum()-1)<1e-8 and fit.x[4]>=fit.x[3]:
            candidates.append(fit)
    if not candidates:
        raise RuntimeError('No feasible successful timing candidate')
    best = min(candidates, key=lambda x: x.fun)
    plan = plateau(); c, hold = plan['concentration'], plan['maintenance']
    expanded = hold*c/(1+c) + np.log1p(c)-np.log1p(c*np.exp(-plan['coast']))
    rows = [dict(regimen='equal',doses=[1/3]*3,times=[0,1,2],J=objective([1/3]*3,[0,1,2])),
            dict(regimen='amounts',doses=fixed.x.tolist(),times=[0,1,2],J=float(-fixed.fun)),
            dict(regimen='timing',doses=best.x[:3].tolist(),times=[0,*best.x[3:].tolist()],J=float(-best.fun)),
            dict(regimen='expanded',doses=[c],times=[0],J=float(expanded),**plan)]
    for i, row in enumerate(rows):
        row['adjacent_gain_percent'] = 0. if i==0 else 100*(row['J']/rows[i-1]['J']-1)
        row['C_upper_bound'] = 1.
        row['injury_upper_bound_full_tail'] = .5
    return dict(seed=SEED,timing_starts=24,successful_timing_starts=len(candidates),rows=rows,
                evidence='Timing is a numerical candidate; expanded design has a supporting analytic certificate for the declared ideal class.')

def trajectory(row, t):
    t=np.asarray(t)
    c=sum(d*np.exp(-np.maximum(t-time,0))*(t>=time) for d,time in zip(row['doses'],row['times']))
    if row['regimen']=='expanded':
        hold=row['maintenance']; plateau_c=row['concentration']
        c=plateau_c*np.where(t<=hold,1.,np.exp(-(t-hold)))
    return c

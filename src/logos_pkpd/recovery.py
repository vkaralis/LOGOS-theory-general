"""Recovery-gated constant-family bounds, parameter stress and noisy information."""
from pathlib import Path
import csv, json
import numpy as np
from scipy.optimize import minimize_scalar, brentq
from scipy.stats import qmc
from scipy.special import logsumexp
from numpy.polynomial.hermite import hermgauss
from . import recovery_models as models

def run(output):
    OUT=Path(output);OUT.mkdir(parents=True,exist_ok=True)
    P=models.P.copy();T0=1.679903;F0=.2
    checks=[]
    def check(name,error,tol):checks.append(dict(check=name,error=float(error),tolerance=tol,passed=bool(error<=tol)))
    def save(name,data): (OUT/(name+'.json')).write_text(json.dumps(data,indent=2,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()))
    def cj(c,p):
     c=np.asarray(c);E=c/(1+c)
     if p['control']=='R_clamped':R=np.ones_like(c)
     else:
      W=p['r']/(1+(c/p['chi'])**p['n']);W=W*0 if p['control']=='recovery_disabled' else W;R=np.divide(W,W+p['alpha']*E,out=np.zeros_like(c,dtype=float),where=(W+p['alpha']*E)>0)
     return p['eta']*E*R-p['rho']*(1-E)
    def cz(c,p):return p['beta']/p['pi']*c/(1+c)*(1+(c/p['chiZ'])**2)
    def bestconstant(p,cap=.22):
     cm=1. if cz(1.,p)<=cap else brentq(lambda c:cz(c,p)-cap,0,1)
     grid=np.linspace(0,cm,10001);vals=cj(grid,p);j=int(np.argmax(vals));lo=grid[max(0,j-1)];hi=grid[min(len(grid)-1,j+1)]
     opt=minimize_scalar(lambda c:-float(cj(c,p)),bounds=(lo,hi),method='bounded',options={'xatol':1e-14});candidates=[(opt.x,-opt.fun),(0.,float(cj(0.,p))),(cm,float(cj(cm,p)))];c,val=max(candidates,key=lambda a:a[1]);return dict(input=c,J=val,injury=float(cz(c,p)),max_admissible_constant=cm)
    constant=bestconstant(P);cycle,_=models.periodic(T0,F0,P)
    # A conservative global upper bound uses |J'| <= L on [0,1]. With n=4,
    # D=r(1+c)+alpha*c*(1+(c/chi)^4); D>=r and |D-cD'|<=D+cD'.
    L=P['eta']*(3+P['alpha']/P['r']*(2+6/P['chi']**4))+P['rho'];grid=np.linspace(0,constant['max_admissible_constant'],1000001)
    upper=float(cj(grid,P).max()+L*(grid[1]-grid[0])/2);lower=cycle['J']-upper
    check('constant local maximum below global derivative-grid upper bound',max(0,constant['J']-upper),1e-12)
    check('cycle positive strict margin against full constant family',max(0,1e-5-lower),0)
    check('cycle satisfies injury cap',max(0,cycle['peak_Z']-.22),1e-9)
    check('cycle periodic closure',cycle['closure_residual'],1e-8)
    nulls=[]
    for label in ['R_clamped','recovery_disabled']:
     p=dict(P,control=label);r=bestconstant(p);c,_=models.periodic(T0,F0,p)
     nulls.append(dict(control=label,constant=r,cycle=c,cycle_minus_constant=c['J']-r['J']))
    save('A_constant_family',dict(comparison='Common average input cap 1, peak-input cap 5, injury cap .22; realized amounts may differ',constant=constant,constant_global_upper_bound=upper,derivative_absolute_bound=L,grid_intervals=1000000,cycle=cycle,certified_numeric_lower_margin=lower,null_controls=nulls))
    print('A full constant-family comparison complete; strict lower margin',lower,flush=True)
    # Deterministic parameter box: not a confidence or probability distribution.
    names=['kappa','alpha','pi','beta','eta','rho','chi','chiZ','r'];factors=.8+.4*qmc.Sobol(len(names),scramble=False).random_base2(5)
    tt,ff=np.meshgrid(np.geomspace(.25,8,9),np.array([.2,.25,.33,.5,.67,.8,1.]));tt=tt.ravel();ff=ff.ravel()
    rows=[]
    for i,fac in enumerate(factors):
     p=dict(P,**{k:P[k]*v for k,v in zip(names,fac)});ref=bestconstant(p)
     j,z=models.screen(tt,ff,p,N=512);eligible=z<=.22-1e-5
     if not eligible.any():selected=None
     else:
      ix=int(np.argmax(np.where(eligible,j,-np.inf)));selected,_=models.periodic(float(tt[ix]),float(ff[ix]),p)
     fixed,_=models.periodic(T0,F0,p)
     row=dict(case=i,parameters={k:p[k] for k in names},constant_J=ref['J'],fixed_J=fixed['J'],fixed_peak_Z=fixed['peak_Z'],fixed_feasible=fixed['peak_Z']<=.22,fixed_increment=fixed['J']-ref['J'],selected_T=selected['T'] if selected else None,selected_f=selected['f'] if selected else None,selected_J=selected['J'] if selected else None,selected_peak_Z=selected['peak_Z'] if selected else None,selected_increment=selected['J']-ref['J'] if selected else None,selected_feasible=bool(selected and selected['peak_Z']<=.22))
     if selected:
      check('screen vs independent selected ODE case '+str(i),abs(j[ix]-selected['J']),8e-6)
      check('selected full-path injury case '+str(i),max(0,selected['peak_Z']-.22),1e-9)
     rows.append(row);print('A robustness',i+1,'/32',flush=True)
    with (OUT/'A_robustness.csv').open('w',newline='') as f:
     writer=csv.DictWriter(f,fieldnames=[k for k in rows[0] if k!='parameters']);writer.writeheader();writer.writerows([{k:v for k,v in r.items() if k!='parameters'} for r in rows])
    save('A_robustness',dict(parameter_box_percent=20,method='32 unscrambled Sobol vectors; deterministic stress test, no probability interpretation',finite_rectangle_grid=63,rows=rows,summary=dict(fixed_feasible=sum(r['fixed_feasible'] for r in rows),fixed_feasible_strict=sum(r['fixed_feasible'] and r['fixed_increment']>1e-6 for r in rows),selected_feasible=sum(r['selected_feasible'] for r in rows),selected_strict=sum(r['selected_feasible'] and r['selected_increment']>1e-6 for r in rows),minimum_selected_increment=min(r['selected_increment'] for r in rows if r['selected_increment'] is not None))))
    # Imperfect information: Bayes-optimal waiting under a fixed three-rate prior.
    rates=np.array([.5,1.,2.]);A=.8;LD=.15;g=.4
    utility=lambda w:A*(-np.expm1(-rates*w))-LD*np.expm1(g*w)
    cal=minimize_scalar(lambda w:-utility(w).mean(),bounds=(0,6),method='bounded',options={'xatol':1e-12});calendar=-cal.fun
    star=np.maximum(0,np.log(A*rates/(LD*g))/(rates+g));ideal=float(np.mean([utility(w)[i] for i,w in enumerate(star)]))
    def info(sd,n):
     if sd==0:return ideal
     nodes,weights=hermgauss(n);total=0.
     for i,r in enumerate(rates):
      for node,weight in zip(nodes,weights):
       y=np.log(r)+np.sqrt(2)*sd*node;logp=-.5*((y-np.log(rates))/sd)**2;p=np.exp(logp-logsumexp(logp))
       opt=minimize_scalar(lambda w:-p@utility(w),bounds=(0,6),method='bounded',options={'xatol':1e-11});total+=weight/np.sqrt(np.pi)/3*utility(opt.x)[i]
     return float(total)
    noise=[]
    for sd in [0,.25,.5,1.,2.]:
     v=info(sd,80);v2=info(sd,160);check('information quadrature convergence SD '+str(sd),abs(v-v2),2e-6)
     noise.append(dict(log_rate_measurement_SD=sd,expected_value=v2,calendar_value=calendar,information_increment=v2-calendar,relative_gain_percent=100*(v2/calendar-1),break_even_cost=v2-calendar,net_gain_at_cost_002=v2-calendar-.002,net_gain_at_cost_01=v2-calendar-.01))
    save('A_imperfect_information',dict(calendar_wait=cal.x,calendar_value=calendar,ideal_value=ideal,prior_rates=rates.tolist(),prior_probabilities=[1/3]*3,shared_wait_bounds=[0,6],noise_results=noise,interpretation='Decision-model value only; sensor errors and costs are hypothetical, not estimated from clinical observations'))
    print('A imperfect information complete',flush=True)
    record=dict(all_passed=all(c['passed'] for c in checks),checks=checks)
    save('A_mathematical_verification',record)
    if not record['all_passed']:raise RuntimeError('Recovery extension verification failed')
    return record

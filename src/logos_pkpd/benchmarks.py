"""Independent analytic, quadrature and ODE benchmark verification."""
from pathlib import Path
import json, sys
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize, minimize_scalar

def run(output):
    OUT=Path(output); OUT.mkdir(parents=True,exist_ok=True)
    checks=[]; results={}
    def check(name,error,tol):
        checks.append(dict(name=name,error=float(error),tolerance=float(tol),passed=bool(error<=tol)))
    
    def emax(doses,times,T=3):
        C=0.;value=0.;prev=0.
        for t,d in zip(times,doses):
            value+=np.log((1+C)/(1+C*np.exp(-(t-prev))))
            C=C*np.exp(-(t-prev))+d;prev=t
        return value+np.log((1+C)/(1+C*np.exp(-(T-prev))))
    j0=emax([1/3]*3,[0,1,2])
    r1=minimize(lambda x:-emax(x,[0,1,2]),[.6,.4,0],method='SLSQP',bounds=[(0,1)]*3,
                constraints={'type':'eq','fun':lambda x:sum(x)-1},options={'ftol':1e-14,'maxiter':500})
    def timing(x):return emax(x[:3],[0,x[3],x[4]])
    r2=minimize(lambda x:-timing(x),[.53,.235,.235,.58,1.17],method='SLSQP',bounds=[(0,1)]*3+[(0,3)]*2,
                constraints=[{'type':'eq','fun':lambda x:sum(x[:3])-1},{'type':'ineq','fun':lambda x:x[4]-x[3]}],options={'ftol':1e-14,'maxiter':500})
    c=brentq(lambda c:c*(4-np.log(2+1/c))-1,.25,1,xtol=1e-14)
    delta=np.log(2+1/c);L=3-delta;j3=L*c/(1+c)+np.log((1+c)/(1+c*np.exp(-delta)))
    vals=[j0,-r1.fun,-r2.fun,j3]
    results['finite_emax']=dict(objectives=vals,adjacent_gain_percent=[100*(b/a-1) for a,b in zip(vals[:-1],vals[1:])],amounts=r1.x.tolist(),timing=r2.x.tolist(),plateau=c,maintenance=L,coast=delta)
    for i,(val,reported) in enumerate(zip(vals,[.629606,.672945,.677043,.679553])):check('Finite reference row '+str(i),abs(val-reported),5e-7)
    
    # Emax marginal support checks at multiple horizons and physical parameters.
    worst=0.;cases=0
    for k,V,E,D in [(1,1,1,1),(.4,2,.7,3),(2,.5,4,.8),(1.7,3,2,4)]:
        threshold=np.log(2+E*V/D)/k
        for T in [threshold*.65,threshold,threshold*1.5,threshold*3]:
            if T<=threshold:
                c0=D/V;hold=0.;C=lambda t:c0*np.exp(-k*t)
            else:
                c0=brentq(lambda z:V*z*(1+k*T-np.log(2+E/z))-D,D/(V*(1+k*T)),D/V)
                hold=T-np.log(2+E/c0)/k;C=lambda t:c0 if t<=hold else c0*np.exp(-k*(t-hold))
            def marg(s):return quad(lambda t:E/(E+C(t))**2*np.exp(-k*(t-s))/V,s,T,points=[hold] if s<hold<T else None,epsabs=1e-11)[0]
            mu=marg(0);worst=max(worst,max(marg(s)-mu for s in np.linspace(0,T,61)));cases+=1
    check('Emax supporting marginal inequality '+str(cases)+' cases',max(0,worst),1e-9)
    
    # Theory B full three-regime condition, tested over arbitrary input times.
    worst=0.;cases=0
    for gamma in [.05,.25,.7,1,1.3]:
        for H in [.05,.2,.7,1,2,4,8]:
            a=4.;b=k=1.;q=gamma*a;T=H;tau=None
            if gamma>=1 or H<=-np.log(1-gamma):
                C=lambda t:0.;support=[]
            else:
                hc=np.arccosh(1/(1-gamma))
                if H<=hc:
                    d=2*(a*(1-np.exp(-H))-q)/(1-np.exp(-2*H));C=lambda t:d*np.exp(-t);support=[0.]
                else:
                    c0=a-q;tau=T-hc;C=lambda t:c0 if t<=tau else c0*np.exp(-(t-tau));support=[0,tau/2,tau]
            def marg(s):return quad(lambda t:np.exp(-(t-s))*(a-C(t)),s,T,points=[tau] if tau is not None and s<tau<T else None,epsabs=1e-11,epsrel=1e-11)[0]
            err=max(0,max(marg(s)-q for s in np.linspace(0,T,101)))
            if support:err=max(err,max(abs(marg(s)-q) for s in support))
            worst=max(worst,err);cases+=1
    check('Theory B supporting condition '+str(cases)+' cases',worst,2e-8)
    a=4.;B=3*(5-np.arccosh(4/3));T=4.;coast=np.arccosh(4/3);tau=T-coast
    J=quad(lambda t:4*(3 if t<=tau else 3*np.exp(-(t-tau)))-.5*(3 if t<=tau else 3*np.exp(-(t-tau)))**2,0,T,points=[tau])[0]
    ref=quad(lambda t:4*B/T*(1-np.exp(-t))-.5*(B/T*(1-np.exp(-t)))**2,0,T)[0]
    eps=.5
    lc=brentq(lambda l:(eps+(4-eps)/np.cosh(l))*(5-l)+eps*l-B,0,4)
    cc=eps+(4-eps)/np.cosh(lc);tc=T-lc
    floor=quad(lambda t:4*(cc if t<=tc else eps+(cc-eps)*np.exp(-(t-tc)))-.5*(cc if t<=tc else eps+(cc-eps)*np.exp(-(t-tc)))**2,0,T,points=[tc])[0]
    results['theory_B_allocation']=dict(B=B,J=J,reference_J=ref,floor_J=floor,shape_gain=floor-ref,stop_gain=J-floor,total_gain=J-ref,coast=coast)
    for name,val,reported in [('B dose',B,12.6139036163),('B optimum',J,28.826262),('B floor',floor,28.686784),('B stopping increment',J-floor,.139478)]:check(name,abs(val-reported),6e-7)
    
    # Closed-form injury models are independently integrated, including equal rates.
    worst=0.;cases=0
    for k in [.3,1,2]:
        for rho in [0,.2,2*k,3*k]:
            for z0 in [0,.3]:
                A0=.7;tt=np.linspace(0,12/k,201)
                if rho==2*k:z=(z0+A0*tt)*np.exp(-2*k*tt)
                else:z=(z0+A0/(2*k-rho))*np.exp(-rho*tt)-A0/(2*k-rho)*np.exp(-2*k*tt)
                sol=solve_ivp(lambda t,y:[A0*np.exp(-2*k*t)-rho*y[0]],(0,tt[-1]),[z0],method='DOP853',t_eval=tt,rtol=1e-11,atol=1e-12)
                worst=max(worst,float(max(abs(sol.y[0]-z))));cases+=1
    check('Quadratic injury tails '+str(cases)+' cases',worst,2e-9)
    tf=np.log(10)/1.8;tp=np.log(5)/.8
    results['delayed_injury']=dict(quadratic_peak_time=tf,quadratic_peak=(np.exp(-.2*tf)-np.exp(-2*tf))/1.8,theory_B_peak_time=tp,theory_B_peak=(np.exp(-.2*tp)-np.exp(-tp))/.8)
    
    # Conditional restart and ideal information benchmark.
    def v(w,r):return .8*(-np.expm1(-r*w))-.15*np.expm1(.4*w)
    wait=np.log(.8/.06)/1.4
    calendar=minimize_scalar(lambda w:-np.mean([v(w,r) for r in [.5,1,2]]),bounds=(0,10),method='bounded',options={'xatol':1e-12})
    ideal=np.mean([v(np.log(.8*r/.06)/(r+.4),r) for r in [.5,1,2]])
    results['theory_A_wait']=dict(wait=wait,value=v(wait,1),calendar_wait=calendar.x,calendar_value=-calendar.fun,ideal_value=ideal,gain_percent=100*(ideal/(-calendar.fun)-1),positive_window_end=brentq(lambda w:v(w,1),1,10))
    check('Theory A wait',abs(wait-1.850191),5e-7);check('Theory A waiting value',abs(v(wait,1)-.509820),5e-7)
    
    payload=dict(purpose='Independent benchmark implementation checks',python=sys.version,numpy=np.__version__,scipy=scipy.__version__,checks=checks,results=results,all_passed=all(c['passed'] for c in checks))
    (OUT/'independent_benchmarks.json').write_text(json.dumps(payload,indent=2)+'\n')
    if not payload['all_passed']:
        raise RuntimeError('Benchmark verification failed')
    return payload

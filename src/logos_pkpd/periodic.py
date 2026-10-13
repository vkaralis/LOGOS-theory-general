"""Dimensionless periodic models evaluated by integrating-factor quadrature."""
import numpy as np
from scipy.integrate import cumulative_simpson, simpson
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar, brentq
def activ(C,sig):return C**3/(1+C**3) if sig else C/(1+C)
def rates(C,sig,dep,r0,rho0):
 q=activ(C,sig);v=C**8/(.35**8+C**8);return q,r0*(1-q)**4 if dep else np.full_like(C,r0),2*q,v,rho0*(1-v) if dep else np.full_like(C,rho0)
def cycle(T,sig=True,dep=True,r0=.7,rho0=1.,N=8193):
 t=np.linspace(0,T,N);C=.5*T/(-np.expm1(-T))*np.exp(-t);q,r,d,v,rho=rates(C,sig,dep,r0,rho0)
 def affine(source,decay):
  integral=cumulative_simpson(decay,x=t,initial=0);forcing=cumulative_simpson(source*np.exp(integral),x=t,initial=0)
  initial=np.exp(-integral[-1])*forcing[-1]/(-np.expm1(-integral[-1]));return np.exp(-integral)*(initial+forcing),initial
 a,a0=affine(r,r+d);z,z0=affine(v,rho)
 spl=CubicSpline(t,z);inds=np.r_[0,np.where((z[1:-1]>=z[:-2])&(z[1:-1]>=z[2:]))[0]+1,N-1];peak=max(z[inds])
 for i in inds:
  if 0<i<N-1:peak=max(peak,-minimize_scalar(lambda tt:-float(spl(tt)),bounds=(t[i-1],t[i+1]),method='bounded',options={'xatol':1e-12}).fun)
 return dict(T=T,J=float(simpson(a*q,x=t)/T),peak_C=float(C[0]),peak_Z=float(peak),mean_Z=float(simpson(z,x=t)/T),closure=max(abs(a[-1]-a0),abs(z[-1]-z0)))
def constant(sig=True,dep=True,r0=.7,rho0=1.):
 q,r,d,v,rho=rates(np.array([.5]),sig,dep,r0,rho0);return dict(J=float(q[0]*r[0]/(r[0]+d[0])),Z=float(v[0]/rho[0]))
def search(sig,dep,r0=.7,rho0=1.):
 ts=np.linspace(.08,12,121);rows=[cycle(float(t),sig,dep,r0,rho0) for t in ts]
 def choose(cap):
  safe=[v for v in rows if v['peak_C']<=4 and v['peak_Z']<=cap];extra=[]
  for i in range(len(ts)-1):
   for constraint in ['peak_C','peak_Z']:
    bound=4 if constraint=='peak_C' else cap
    if (rows[i][constraint]-bound)*(rows[i+1][constraint]-bound)<0:
     root=brentq(lambda tt:cycle(tt,sig,dep,r0,rho0)[constraint]-bound,ts[i],ts[i+1],xtol=2e-9);extra.append(cycle(root,sig,dep,r0,rho0))
  safe += [v for v in extra if v['peak_C']<=4+1e-8 and v['peak_Z']<=cap+1e-8]
  if not safe:return None
  best=max(safe,key=lambda v:v['J']);i=int(np.argmin(abs(ts-best['T'])))
  opt=minimize_scalar(lambda tt:-cycle(tt,sig,dep,r0,rho0)['J'],bounds=(ts[max(0,i-1)],ts[min(len(ts)-1,i+1)]),method='bounded',options={'xatol':2e-9})
  v=cycle(opt.x,sig,dep,r0,rho0)
  if v['peak_C']<=4+1e-8 and v['peak_Z']<=cap+1e-8 and v['J']>best['J']:best=v
  return best
 return choose,rows

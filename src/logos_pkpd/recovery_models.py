"""Dimensionless recovery-gated rectangular input model."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.linalg import expm
P=dict(kappa=4.,alpha=2.,pi=.7,beta=.2,eta=1.,rho=.08,chi=.35,chiZ=1.,n=4.,r=1.,gate='hill',control='full')

def coeff(C,p):
    E=C/(1+C)
    if p['gate']=='hill': W=1/(1+(C/p['chi'])**p['n'])
    else: W=np.exp(-np.log(2)*(C/p['chi'])**(p['n']/(2*np.log(2))))
    sr=p['r']*W
    if p['control']=='recovery_disabled':sr=np.zeros_like(C)
    dr=sr+p['alpha']*E
    wz=np.ones_like(C) if p['control']=='ordinary_repair' else 1/(1+(C/p['chiZ'])**2)
    return E,sr,dr,p['beta']*E,p['pi']*wz

def screen(T,f,p=P,N=256):
    T,f=np.broadcast_arrays(np.atleast_1d(T).astype(float),np.atleast_1d(f).astype(float))
    amp=1/f;on=T*f;off=T-on;k=p['kappa']
    C0=amp*(-np.expm1(-k*on))*np.exp(-k*off)/(-np.expm1(-k*T))
    Cswitch=amp+(C0-amp)*np.exp(-k*on)
    AR=np.ones_like(T);BR=np.zeros_like(T);AZ=np.ones_like(T);BZ=np.zeros_like(T)
    segments=[]
    for duration,initial,u in [(on,C0,amp),(off,Cswitch,0)]:
        dt=duration/N
        uvec=np.broadcast_to(u,T.shape)
        c=uvec[:,None]+(initial-uvec)[:,None]*np.exp(-k*dt[:,None]*(np.arange(N)+.5))
        E,sr,dr,sz,dz=coeff(c,p)
        rr=np.exp(-dr*dt[:,None]);rz=np.exp(-dz*dt[:,None])
        br=sr*(-np.expm1(-dr*dt[:,None]))/dr
        bz=sz*(-np.expm1(-dz*dt[:,None]))/dz
        segments.append((dt,E,sr,dr,sz,dz,rr,rz,br,bz))
        for j in range(N):
            AR*=rr[:,j];BR=rr[:,j]*BR+br[:,j]
            AZ*=rz[:,j];BZ=rz[:,j]*BZ+bz[:,j]
    R=BR/(1-AR);Z=BZ/(1-AZ);peak=Z.copy();integral=np.zeros_like(T)
    for dt,E,sr,dr,sz,dz,rr,rz,br,bz in segments:
        for j in range(N):
            mid=np.exp(-dr[:,j]*dt/2)*R+sr[:,j]*(-np.expm1(-dr[:,j]*dt/2))/dr[:,j]
            if p['control']=='R_clamped':mid=np.ones_like(mid)
            integral+=dt*(p['eta']*E[:,j]*mid-p['rho']*(1-E[:,j]))
            R=rr[:,j]*R+br[:,j];Z=rz[:,j]*Z+bz[:,j];peak=np.maximum(peak,Z)
    return integral/T,peak

def reference(p=P):
    E,sr,dr,sz,dz=coeff(np.array([1.]),p)
    R=sr/dr
    if p['control']=='R_clamped':R[:]=1
    return float(p['eta']*E[0]*R[0]-p['rho']*(1-E[0])),float(sz[0]/dz[0])

def periodic(T,f,p=P,method='DOP853',rtol=2e-11,absorption=False):
    on=T*f;u=1/f;k=p['kappa'];segments=[(0,on,u),(on,T,0.)]
    if absorption:
        ka=12.;M=np.array([[-ka,0],[k,-k]]);B=np.array([ka,0.])
        def step(x,v,t):
            A=expm(M*t);return A@x+np.linalg.solve(M,(A-np.eye(2))@B*v)
        A=expm(M*T);b=step(step(np.zeros(2),u,on),0,T-on)
        x0=np.linalg.solve(np.eye(2)-A,b);xs=step(x0,u,on)
        def concentration(t):
            return float(step(x0,u,t)[1]) if t<=on else float(step(xs,0,t-on)[1])
    else:
        c0=u*(-np.expm1(-k*on))*np.exp(-k*(T-on))/(-np.expm1(-k*T))
        cs=u+(c0-u)*np.exp(-k*on)
        def concentration(t):return u+(c0-u)*np.exp(-k*t) if t<=on else cs*np.exp(-k*(t-on))
    # One-cycle affine map integrated separately from the screening algorithm.
    def affine(t,y):
        E,sr,dr,sz,dz=coeff(np.array([concentration(t)]),p)
        return [-dr[0]*y[0],sr[0]-dr[0]*y[1],-dz[0]*y[2],sz[0]-dz[0]*y[3]]
    y=np.array([1.,0,1,0])
    for start,end,_ in segments:
        so=solve_ivp(affine,(start,end),y,method=method,rtol=rtol,atol=rtol/100,max_step=T/100)
        assert so.success;y=so.y[:,-1]
    R0=y[1]/(1-y[0]);Z0=y[3]/(1-y[2])
    def fun(t,y):
        E,sr,dr,sz,dz=coeff(np.array([concentration(t)]),p)
        R=1 if p['control']=='R_clamped' else y[0]
        return [sr[0]-dr[0]*y[0],sz[0]-dz[0]*y[1],p['eta']*E[0]*R-p['rho']*(1-E[0])]
    y=np.array([R0,Z0,0.]);sols=[];extrema=[]
    for start,end,_ in segments:
        so=solve_ivp(fun,(start,end),y,method=method,rtol=rtol,atol=rtol/100,max_step=T/150,dense_output=True)
        assert so.success
        points=np.linspace(start,end,1501);der=np.array([fun(t,so.sol(t))[1] for t in points])
        roots=[brentq(lambda t:fun(t,so.sol(t))[1],a,b,xtol=1e-13) for a,b,da,db in zip(points[:-1],points[1:],der[:-1],der[1:]) if da*db<0]
        extrema.extend(float(so.sol(t)[1]) for t in [start,end]+roots)
        y=so.y[:,-1];sols.append(so)
    ts=np.linspace(0,T,1201);ys=np.array([sols[0 if t<=on else 1].sol(t) for t in ts])
    C=np.array([concentration(t) for t in ts])
    trace=dict(time=ts,u=np.where(ts<on,u,0.),C=C,R=ys[:,0],Z=ys[:,1],D_relative=np.exp(-ys[:,2]))
    return dict(T=float(T),f=float(f),J=float(y[2]/T),peak_Z=max(extrema),
                closure_residual=float(max(abs(y[0]-R0),abs(y[1]-Z0))),mean_C=float(np.trapezoid(C,ts)/T)),trace


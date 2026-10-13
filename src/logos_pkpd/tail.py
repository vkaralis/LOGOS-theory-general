"""Exact cessation injury including zero repair and coincident rates."""
import numpy as np

def injury(t, c0=1., z0=0., alpha=1., k=1., rho=.2):
    t=np.asarray(t,dtype=float)
    if np.any(t<0) or min(c0,z0,alpha,rho)<0 or k<=0:
        raise ValueError('Require t,c0,z0,alpha,rho >= 0 and k > 0')
    drive=alpha*c0*c0
    if np.isclose(rho,2*k,rtol=0,atol=1e-12):
        return (z0+drive*t)*np.exp(-rho*t)
    if rho > 2*k:
        forcing = drive*np.exp(-2*k*t)*(-np.expm1(-(rho-2*k)*t))/(rho-2*k)
    else:
        forcing = drive*np.exp(-rho*t)*(-np.expm1(-(2*k-rho)*t))/(2*k-rho)
    return z0*np.exp(-rho*t)+forcing

def peak(c0=1., z0=0., alpha=1., k=1., rho=.2):
    injury([0],c0,z0,alpha,k,rho)
    drive=alpha*c0*c0
    if rho==0:
        return dict(time=None if drive else 0., value=z0+drive/(2*k), attained=not bool(drive))
    candidates=[(0.,z0)]
    if drive:
        if np.isclose(rho,2*k,rtol=0,atol=1e-12):
            time=1/rho-z0/drive
        else:
            A=drive/(2*k-rho);den=rho*(z0+A)
            ratio=2*k*A/den if den else -1
            time=np.log(ratio)/(2*k-rho) if ratio>0 else -1
        if time>0:
            candidates.append((float(time),float(injury(time,c0,z0,alpha,k,rho))))
    time,value=max(candidates,key=lambda x:x[1])
    return dict(time=time,value=value,attained=True)

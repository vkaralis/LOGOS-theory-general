"""Static Hill-response relaxed concave envelope; not an attainable schedule."""
import numpy as np

def hill(c, h=3., ec50=1., emax=1.):
    c=np.asarray(c,dtype=float)
    if h<=1 or ec50<=0 or emax<=0 or np.any(c<0):
        raise ValueError('Require h>1, ec50/emax>0 and c>=0')
    x=(c/ec50)**h
    return emax*x/(1+x)

def envelope(mean, cap, h=3., ec50=1., emax=1.):
    mean=np.asarray(mean,dtype=float)
    hill(mean,h,ec50,emax)
    if cap<=0 or np.any(mean>cap):raise ValueError('Mean must be inside [0, cap]')
    support=min(cap,ec50*(h-1)**(1/h))
    return np.where(mean<=support,mean/support*hill(support,h,ec50,emax),hill(mean,h,ec50,emax))

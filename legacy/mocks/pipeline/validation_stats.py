"""Predeclared absolute ensemble statistics and covariance-aware shape test."""
import numpy as np
from scipy.stats import t, f

def absolute_statistics(values,target=0.):
    x=np.asarray(values,float); n=len(x)
    mean=float(x.mean()); sem=float(x.std(ddof=1)/np.sqrt(n))
    bound=float(abs(mean-target)+t.ppf(.975,n-1)*sem)
    return dict(n=n,mean=mean,sem=sem,target=float(target),residual=mean-target,bound95=bound)

def slope_statistics(amplitudes,truth=(0,.5,1,2)):
    x=np.asarray(truth,float); y=np.asarray(amplitudes,float)
    xc=x-x.mean(); slopes=y@xc/np.dot(xc,xc)
    result=absolute_statistics(slopes,1.)
    result['slopes']=slopes.tolist(); result['intercepts']=(y.mean(axis=1)-slopes*x.mean()).tolist()
    return result

def hotelling_shape(values):
    x=np.asarray(values,float); n,b=x.shape
    contrasts=np.eye(b)[1:]-np.eye(b)[0]
    y=x@contrasts.T; p=b-1
    cov=np.cov(y,rowvar=False); mu=y.mean(axis=0)
    if n<=p or np.linalg.matrix_rank(cov)<p:
        return dict(p_value=0.,T2=float('nan'),reason='insufficient nonsingular ensemble covariance')
    tsq=float(n*mu@np.linalg.solve(cov,mu))
    return dict(T2=tsq,F=(n-p)*tsq/(p*(n-1)),df=[p,n-p],p_value=float(f.sf((n-p)*tsq/(p*(n-1)),p,n-p)),covariance=cov.tolist(),means=x.mean(axis=0).tolist(),sem=(x.std(axis=0,ddof=1)/np.sqrt(n)).tolist())

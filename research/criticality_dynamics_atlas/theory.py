"""Minimal nulls, all passed as input series to the official MATLAB runner."""
import numpy as np
from scipy.signal import lfilter
from metrics import population_metrics

SEEDS=(3141,3142,3143)
DT=.04
N=24000
BURN=3000

def ar(m,n,noise):
    return lfilter([1],[1,-m],noise)[:n]


def branch(m,drive,rng):
    y=np.zeros(len(drive))
    for i in range(1,len(y)):
        y[i]=rng.poisson(m*y[i-1]+drive[i])
    return y


def markov(n,cut,length,rng):
    # Two-state chain with stationary OFF probability cut and mean OFF length length.
    off=False
    mask=np.ones(n)
    for i,u in enumerate(rng.random(n)):
        off=not (u<1/length) if off else u<cut/(length*(1-cut))
        mask[i]=not off
    return mask


def scenarios():
    rows=[]
    for seed in SEEDS:
        rng=np.random.default_rng(seed)
        eps=rng.normal(size=(3,N+BURN))
        ar1=[(f"ar_{m}",ar(m,N+BURN,eps[0]),dict(m=m,tau=-DT/np.log(m),
                var=1/(1-m*m),chi=1/(1-m))) for m in (.2,.5,.8,.92,.97,.99)]
        variants=list(ar1)
        rows += [(name,x[BURN:],seed,p,1) for name,x,p in ar1]
        # Each constituent has unit stationary variance.
        for slow in (.92,.98,.995):
            for weight in (.05,.25,.65):
                for noise in (0,.40):
                    fast=.45
                    x=(np.sqrt(1-weight)*ar(fast,N+BURN,eps[0]*np.sqrt(1-fast**2))+
                       np.sqrt(weight)*ar(slow,N+BURN,eps[1]*np.sqrt(1-slow**2)))
                    x=np.sqrt(1-noise)*x+np.sqrt(noise)*eps[2]
                    variants.append((f"mix_{slow}_{weight}_{noise}",x,
                        dict(slow=slow,weight=weight,noise=noise,long_tau=-DT/np.log(slow),fast=fast,
                             analytic_tau_int=DT*(.5+(1-noise)*((1-weight)*fast/(1-fast)+weight*slow/(1-slow))))))
        for m in (.2,.5,.8,.92,.97):
            x=branch(m,np.full(N+BURN,10*(1-m)),rng)
            for sampling in (1.,.3,.1):
                y=rng.binomial(x.astype(int),sampling).astype(float)
                variants.append((f"branch_{m}_{sampling}",y,dict(m=m,tau=-DT/np.log(m),sampling=sampling)))
        latent=branch(.8,np.full(N+BURN,2.),rng)
        drive=ar(.97,N+BURN,eps[0]*np.sqrt(1-.97**2))
        for fraction in (.15,.35):
            for length in (3,15):
                y=latent*markov(len(latent),fraction,length,rng)
                variants.append((f"gate_{fraction}_{length}",y,dict(latent_m=.8,off_fraction=fraction,
                                                                    mean_off_s=length*DT)))
        variants.extend([("latent_ungated",latent,dict(latent_m=.8)),
          ("common_poisson",rng.poisson(10*np.exp(.5*drive-.125)).astype(float),dict(latent_m=0.)),
          ("independent_poisson",rng.poisson(10,len(latent)).astype(float),dict(latent_m=0.))])
        rows += [(name,x[BURN:],seed,p,5) for name,x,p in variants]
    return rows


def prepare():
    stats=[]
    fits={1:[],5:[]}
    for name,series,seed,p,order in scenarios():
        i=f"{name}-{seed}-p{order}"
        m=dict(id=i,session=i,region="Null",context=name,state="Synthetic",
               variant="reference" if order==5 else "ar1",seed=seed)
        stats.append(dict(id=i,scenario=name,seed=seed,order=order,truth=p,
                          **population_metrics([series],DT,False)))
        fits[order].append((m,[series]))
    return stats,fits

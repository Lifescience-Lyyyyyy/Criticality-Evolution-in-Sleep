"""Segment-aware dynamics for continuous bins; time is in seconds."""
import numpy as np
from scipy.fft import irfft, next_fast_len, rfft
from scipy.optimize import least_squares
from scipy.signal import find_peaks, welch


def lag_cov(blocks, k, local=True):
    """Common denominator at each lag, pooled only over valid pairs."""
    blocks = [np.asarray(x, float) for x in blocks]
    center = sum(x.sum(0) for x in blocks)/sum(map(len, blocks))
    n = np.zeros(k+1)
    c = np.zeros((k+1,)+blocks[0].shape[1:])
    for x in blocks:
        x = x-(x.mean(0) if local else center)
        d = min(k+1,len(x))
        fft = rfft(x,next_fast_len(2*len(x)-1),axis=0)
        c[:d] += irfft(np.abs(fft)**2,axis=0)[:d]
        n[:d] += len(x)-np.arange(d)
    shape=(k+1,)+(1,)*(c.ndim-1)
    return np.divide(c,n.reshape(shape),out=np.full_like(c,np.nan),where=n.reshape(shape)>0)


def exp_fit(y, dt, start, end, offset):
    """Grid fit with linear amplitude/offset, using inclusive positive lags."""
    a=max(1,round(start/dt)); b=min(len(y),round(end/dt)+1)
    t=np.arange(a,b)*dt
    y=np.asarray(y[a:b],float)
    if len(t)<4 or not np.isfinite(y).all():
        return {}
    tau=np.geomspace(dt/2,max(t)*4,90)
    x=np.exp(-t[None,:]/tau[:,None])
    if offset:
        u=x-x.mean(1,keepdims=True)
        v=y-y.mean()
    else:
        u=x
        v=y
    amp=(u@v)/(u*u).sum(1)
    c=y.mean()-amp*x.mean(1) if offset else np.zeros(len(tau))
    se=np.square(y-amp[:,None]*x-c[:,None]).sum(1)
    j=np.argmin(np.where(amp>0,se,np.inf))
    r2=1-se[j]/max(np.square(y-y.mean()).sum(),1e-15)
    ans=dict(tau=float(tau[j]),amp=float(amp[j]),offset=float(c[j]),
             r2=float(r2),lo=float(t[0]),hi=float(t[-1]),rmse=float(np.sqrt(se[j]/len(t))))
    # Floor/ceiling and poor-decay estimates are retained as diagnostics only.
    ans["valid"]=bool(0<j<len(tau)-1 and amp[j]>.025 and r2>=.35 and tau[j]<=t[-1]*2)
    return ans


def integrated(rho,dt):
    """One-sided positive-lobe integral, with a 2 s maximum horizon."""
    p=np.asarray(rho[1:round(2/dt)+1])
    first=np.flatnonzero(~np.isfinite(p)|(p<=0))
    n=first[0] if len(first) else len(p)
    return float(dt*(.5+p[:n].sum())),not len(first)


def mr_coeff(blocks,k):
    """Pooled regression of x[t+k] on x[t], with per-block centering."""
    out=np.ones(k+1)
    for lag in range(1,k+1):
        cov=var=0.
        for x in blocks:
            if len(x)<=lag:
                continue
            a,b=x[:-lag],x[lag:]
            a=a-a.mean()
            cov+=float(a@(b-b.mean()))
            var+=float(a@a)
        out[lag]=cov/var if var>0 else np.nan
    return out


def mr(blocks,dt):
    c=mr_coeff(blocks,max(8,round(2/dt)))
    full=exp_fit(c,dt,dt,2,False)
    short=exp_fit(c,dt,dt,max(.4,4*dt),False)
    for f in (full,short):
        if f:
            f["m"]=float(np.exp(-dt/f["tau"]))
    return full,short,c


def moments(blocks):
    count=sum(map(len,blocks))
    mean=sum(x.sum(0) for x in blocks)/count
    local=[x-x.mean(0) for x in blocks]
    cov=sum(x.T@x for x in local)/count
    return mean, cov


def pr(c):
    power=np.square(c).sum()
    return float(c.trace()**2/power) if power>0 else None


def unit_metrics(blocks,dt,seed,detail):
    means,c=moments(blocks)
    n=len(means)
    active=c.diagonal()>0
    z=np.sqrt(c.diagonal()[active])
    raw=pr(c)
    cc=c[np.ix_(active,active)]/np.outer(z,z)
    population=sum(x.sum(0) for x in c)
    eig=np.linalg.eigvalsh(c)[::-1].clip(0)
    out=dict(pr_raw=raw,pr_raw_n=raw/n if raw is not None else None,
             pr_corr=pr(cc),pr_corr_n=pr(cc)/len(z) if len(z) else None,
             active_units=int(active.sum()),top_fraction=float(eig[0]/eig.sum()) if eig.sum() else None,
             pair_r=float(cc[np.triu_indices(len(z),1)].mean()) if len(z)>1 else None,
             rate_cv=float(means.std()/means.mean()) if means.mean() else None)
    var_other=population+c.diagonal()-2*c.sum(0)
    denom=np.sqrt(np.maximum(0,c.diagonal()*var_other))
    usable=denom>0
    out["coupling"]=float(np.median((c.sum(0)-c.diagonal())[usable]/denom[usable])) if usable.any() else None
    if not detail:
        return out
    k=max(5,round(8/dt))
    ac=lag_cov(blocks,k)
    pop=lag_cov([x.sum(1) for x in blocks],k)
    auto=ac.sum(1)
    cross=pop-auto
    for lag in [0,1,min(k,round(.4/dt))]:
        pieces=[(x-x.mean(0)) for x in blocks if len(x)>lag]
        direct=sum((x[:len(x)-lag].T@x[lag:]) for x in pieces)/sum(len(x)-lag for x in pieces)
        assert np.allclose([direct.trace(),direct.sum()],[auto[lag],pop[lag]],rtol=2e-7,atol=1e-9)
    rho=np.divide(ac,ac[0],out=np.full_like(ac,np.nan),where=ac[0]>0)
    neurons=[]
    for i in range(n):
        f=exp_fit(rho[:,i],dt,2*dt,.8,True) if means[i]*sum(map(len,blocks))>=100 else {}
        tau_i,tau_int=None,None
        if f.get("valid"):
            tau_i=f["tau"]
        if f and np.isfinite(rho[0,i]):
            tau_int=integrated(rho[:,i],dt)[0]
        neurons.append([float(means[i]/dt),float(c[i,i]/means[i]) if means[i]>0 else None,tau_i,tau_int])
    valid=[u[2] for u in neurons if u[2] is not None]
    rng=np.random.default_rng(seed)
    fixed={}
    for size in (5,10,20,40):
        if size>n:
            continue
        samples=[rng.choice(n,size,False) for _ in range(12)]
        for s in (size,):
            vals=[]
            for idx in samples:
                part=c[np.ix_(idx,idx)]
                if pr(part) is not None:
                    vals.append(pr(part)/s)
            if vals:
                fixed[str(s)]=[float(np.quantile(vals,v)) for v in (.25,.5,.75)]
    # Use a fixed 2 s horizon. When cancellation is severe report the raw contributions.
    h=round(2/dt)
    total=float(np.sum(pop[1:h+1])*dt)
    out.update(neurons=neurons,neuron_tau_median=float(np.median(valid)) if valid else None,
               neuron_tau_iqr=np.quantile(valid,[.25,.75]).tolist() if valid else None,
               neuron_tau_fraction=len(valid)/n,
               neuron_tau_over_250ms=float(np.mean(np.asarray(valid)>.25)) if valid else None,
               neuron_fano_median=float(np.median([u[1] for u in neurons if u[1] is not None])),
               fixed_n=fixed,eigenspectrum=(eig/eig.sum()).tolist(),
               cross_fraction0=float(cross[0]/pop[0]),
               cross_integral=float(np.sum(cross[1:h+1])*dt),
               cross_fraction_integrated=float(np.sum(cross[1:h+1])*dt/total) if abs(total)>.01*pop[0] else None,
               cov_auto=auto.tolist(),cov_cross=cross.tolist(),cov_total=pop.tolist())
    if len(blocks)>1:
        _,s=moments(blocks[:1])
        out["pr_first_window_n"]=pr(s)/n if pr(s) is not None else None
    if n>=10:
        idx=rng.choice(n,n//2,False)
        half=[x[:,idx].sum(1) for x in blocks]
        f,_,_=mr(half,dt)
        out["m_half"]=f.get("m") if f.get("valid") else None
    return out


def lengths(mask):
    e=np.flatnonzero(np.diff(np.r_[False,mask,False]))
    return (e[1::2]-e[::2])[(e[::2]>0)&(e[1::2]<len(mask))]


def events(blocks,dt,kind,cut=0):
    masks=[(x<=cut) if kind=="off" else (x>cut) for x in blocks]
    sizes=np.concatenate([lengths(v) for v in masks])*dt
    spans=np.concatenate([lengths(~v) for v in masks])*dt
    interval=[]
    for v in masks:
        starts=1+np.flatnonzero(np.diff(v.astype(np.int8))==1)
        interval.extend(np.diff(starts)*dt)
    s=np.array([.04,.08,.12,.20,.32,.50,.80,1.28,2.,4.,8.,16.,np.inf])
    return dict(fraction=float(sum(v.sum() for v in masks)/sum(map(len,masks))),
                events=int(len(sizes)),event_hz=float(len(sizes)/(sum(map(len,masks))*dt)),
                transition_hz=float(sum(np.count_nonzero(np.diff(v)) for v in masks)/(sum(map(len,masks))*dt)),
                mean_s=float(sizes.mean()) if len(sizes) else None,
                p95_s=float(np.quantile(sizes,.95)) if len(sizes) else None,
                cv=float(sizes.std()/sizes.mean()) if len(sizes)>1 else None,
                active_mean_s=float(spans.mean()) if len(spans) else None,
                interevent_s=float(np.mean(interval)) if interval else None,
                hist=np.histogram(sizes,np.r_[0,s])[0].tolist())


def spectrum(blocks,dt):
    n=round(8/dt)
    xx=[x for x in blocks if len(x)>=n]
    if not xx:
        return {}
    sp=[welch(x,fs=1/dt,nperseg=n,noverlap=n//2) for x in xx]
    f=sp[0][0]
    p=np.average([q[1] for q in sp],axis=0,weights=[1+(len(x)-n)//(n-n//2) for x in xx])
    p_ok=(p>0)&(f>=.375)&(f<=min(10,.4/dt))
    x,y=f[p_ok],np.log(p[p_ok])
    knee=slope=None
    if len(x)>10:
        slope=float(np.polyfit(np.log(x),y,1)[0])
        approx=np.polyval(np.polyfit(np.log(x),y,1),np.log(x))
        def err(v):
            return v[0]-np.log1p((x/np.exp(v[1]))**v[2])-y
        fit=least_squares(err,[max(y),0,1],bounds=([-60,np.log(.25),.1],[60,np.log(8),5]),
                          loss="soft_l1",f_scale=.25,max_nfev=60)
        loc=float(np.exp(fit.x[1]))
        rss=np.square(err(fit.x)).sum()
        alt=np.square(approx-y).sum()
        if f[1]<loc<x.max()*.8 and loc>x.min()*1.2 and rss>0 and alt>rss:
            if len(x)*np.log(alt/rss)>6:
                knee=loc
    lo=p[(f>=.25)&(f<=2)].sum()*f[1]
    hi=p[(f>=4)&(f<=10)].sum()*f[1] if f[-1]>=10 else None
    peaks,_=find_peaks(p,prominence=max(p[1:])*.15)
    return dict(freq=f.tolist(),psd=p.tolist(),slope=slope,knee_hz=knee,
                low_power=float(lo),high_power=float(hi) if hi is not None else None,
                low_high=float(lo/hi) if hi else None,peak_hz=f[peaks[:8]].tolist())


def population_metrics(blocks,dt,profiles=True):
    x=np.concatenate(blocks).astype(float)
    if not x.var():
        return {}
    cov=lag_cov(blocks,round(8/dt))
    rho=cov/cov[0]
    tau=exp_fit(rho,dt,dt,2,True)
    f,s,m=mr(blocks,dt)
    tau_i,censored=integrated(rho,dt)
    mean=x.mean()
    silence=events(blocks,dt,"off")
    burst=events(blocks,dt,"high",float(np.quantile(x,.9)))
    out=dict(tau_population=tau["tau"] if tau.get("valid") else None,tau_fit=tau,
             tau_integrated=tau_i,tau_int_truncated=censored,m_MR=f["m"] if f.get("valid") else None,
             tau_MR=f["tau"] if f.get("valid") else None,m_short=s["m"] if s.get("valid") else None,
             mr_fit=f,mean_rate=float(mean/dt),variance=float(x.var()),
             fano=float(x.var()/mean) if mean>0 else None,cv=float(x.std()/mean) if mean>0 else None,
             silent_fraction=silence["fraction"],off_duration=silence["p95_s"],
             off_mean=silence["mean_s"],off_hz=silence["event_hz"],
             block_mean_cv=float(np.std([y.mean() for y in blocks])/abs(mean)) if mean else None,
             lag1=float(rho[1]),persist1s=float(rho[round(1/dt)]),persist4s=float(rho[round(4/dt)]),
             burst_p95=burst["p95_s"])
    if profiles:
        out.update(acf=rho.tolist(),mr_coeff=m.tolist(),off=silence,burst=burst,**spectrum(blocks,dt))
    return out


def ar_modes(coeff,dt):
    lam=np.roots(np.r_[1,-np.asarray(coeff)])
    return [{"tau":float(-dt/np.log(abs(z))),
             "period":float(2*np.pi*dt/abs(np.angle(z))) if np.angle(z)!=0 else None}
            for z in sorted(lam,key=abs,reverse=True) if 0<abs(z)<1]


def shift(x,rng):
    n=len(x)
    return np.column_stack([np.roll(col,int(rng.integers(25,n-24))) for col in x.T])


def iaaft(x,rng,steps=100):
    target=np.abs(rfft(x))
    y=rng.permutation(x)
    ordered=np.sort(x)
    for _ in range(steps):
        phase=np.angle(rfft(y))
        z=irfft(target*np.exp(1j*phase),n=len(x))
        y[np.argsort(z,kind="stable")]=ordered
    return y

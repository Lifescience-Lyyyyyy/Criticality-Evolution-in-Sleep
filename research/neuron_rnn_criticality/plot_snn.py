"""Plot the measured physiology, timescales, shifts and pulse responses."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

LAB={"figure7":"Fig. 7 · synchronous kicks","figure7_supp2":"Fig. 7 suppl. 2 · independent noise"}
CT={"figure7":"#007884","figure7_supp2":"#8C4A1C"}
CELL={"E":"#9E321F","I":"#245C9A"}
PLT={"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"font.size":8,
     "pdf.fonttype":42,"axes.spines.right":False,"axes.spines.top":False,"axes.titleweight":"bold"}
plt.rcParams.update(PLT)


def read(path):
    return json.loads(path.read_text())


def grid(r,c,size):
    f,axes=plt.subplots(r,c,figsize=size,layout="constrained",squeeze=False)
    for i,ax in enumerate(axes.flat):
        ax.text(-.085,1.025,chr(65+i),transform=ax.transAxes,weight="bold")
    return f,list(axes.flat)


def save(fig,out,index,name):
    fig.savefig(out/f"{index:02d}_{name}.pdf")
    plt.close(fig)


def lineband(ax,xx,curves,color,label,linestyle="-"):
    y=np.asarray(curves)
    if not len(y):
        return
    ax.plot(xx,np.median(y,axis=0),color=color,lw=1.4,ls=linestyle,label=label)
    if len(y)>1:
        ax.fill_between(xx,y.min(0),y.max(0),color=color,alpha=.12,lw=0)


def scope(fig,caption):
    fig.supxlabel(caption,fontsize=8)


def plot(data,phys,fit,out,pulses=None):
    fs={x["id"]:x for x in fit["records"]}
    conditions=[k for k in LAB if any(x["manifest"]["condition"]==k for x in data["runs"])]
    refs={k:[z for z in data["runs"] if z["manifest"]["condition"]==k] for k in conditions}
    whole="Original author parameters · Jercog network, N=5,000 (4,000 E + 1,000 I) · three 100 s runs per condition"
    pair=tau_runs=list(data["runs"])
    fig,axes=grid(len(conditions),2,(9.3,2.22*len(conditions)+.6))
    for i,key in enumerate(conditions):
        a,b=axes[2*i:2*i+2]
        q=next(r for r in phys["runs"] if r["manifest"]["condition"]==key)
        m,z=q["manifest"],q["preview"]
        for k,offset in (("E",0),("I",65)):
            t,j=np.asarray(z["raster_"+k]).T
            a.scatter(t,j+offset,s=1.4,c=CELL[k],rasterized=False,linewidths=0)
        t=np.asarray(z["time_s"])
        msk=np.asarray(z["up"],bool)
        edges=np.r_[0,np.flatnonzero(np.diff(msk))+1,len(msk)]
        for j,k in zip(edges[:-1],edges[1:]):
            if msk[j]:
                b.axvspan(t[j],t[k-1]+.01,color="#D8EADC",lw=0,zorder=0)
        for k in ("E","I"):
            b.plot(t,np.convolve(z["rate_"+k],np.ones(5)/5,mode="same"),lw=1.0,c=CELL[k],label=k)
        a.set(title=LAB[key],ylabel="Recorded cell index",ylim=(-2,132),xlim=(12,18),xlabel="Simulation time (s)")
        b.set(title=f"Delay seed {m['connectivity_seed']}; noise seed {m['simulation_seed']}",
              ylabel="Population rate (Hz/cell)",xlim=(12,18),xlabel="Simulation time (s)")
        y=b.get_ylim()[1]
        b.set_ylim(top=y*1.27)
        b.legend(handles=[Line2D([],[],color=CELL[k],label=k) for k in ("E","I")]+
                         [Line2D([],[],lw=5,color="#D8EADC",label="UP")],
                 loc="upper right",ncol=3,frameon=False,fontsize=7)
    scope(fig,whole)
    save(fig,out,1,"original_rasters")

    fig,axes=grid(2,2,(9.3,6.2))
    qgrid=np.geomspace(.035,1.6,90)
    for ax,k in zip(axes[:2],conditions):
        for t,sl in (("E",slice(0,4000)),("I",slice(4000,5000))):
            curves=[]
            fra=[]
            for r in refs[k]:
                x=np.asarray(r["individual"],dtype=float)[sl,1]
                use=np.sort(x[np.isfinite(x)])
                fra.append(len(use)/len(x))
                if len(use):
                    curves.append(np.searchsorted(use,qgrid,side="right")/len(use))
            lineband(ax,qgrid,curves,CELL[t],f"{t}, median fit rate {np.median(fra):.0%}")
        ax.set(xscale="log",ylim=(0,1.24),xlabel="Single-cell exponential τ (s)",
               ylabel="Fraction of identifiable units",title=LAB[k])
        ax.legend(frameon=False,loc="upper left",fontsize=7)
    for k in conditions:
        rows=refs[k]
        lineband(axes[2],np.arange(151)*.04,[x["population"]["acf"][:151] for x in rows],CT[k],LAB[k])
        for r in rows:
            x=np.asarray(r["individual"],dtype=float)[:,1]
            p=r["population"]["tau_population"]
            if p is not None and np.isfinite(x).any():
                axes[3].scatter(np.nanmedian(x),p,s=28,c=CT[k])
    p=axes[3]
    if p.has_data():
        bounds=np.array([p.get_xlim(),p.get_ylim()])
        low=bounds[:,0].min()
        high=bounds[:,1].max()
        p.plot([low,high],[low,high],":",c="#777777")
    p.set(xlabel="Median single-cell τ (s)",ylabel="Population exponential τ (s)")
    axes[2].set(xlabel="Lag (s)",ylabel="Population normalized ACF",ylim=(None,1.22))
    axes[2].legend(frameon=False,fontsize=7)
    scope(fig,whole+" · 90 s measured after burn-in")
    save(fig,out,2,"timescales")

    fig,axes=grid(1,len(conditions),(9.3,3.05))
    palette={"auto":"#2479A1","cross":"#8E4635","total":"#21354B"}
    for ax,k in zip(axes,conditions):
        for component in ("auto","cross","total"):
            y=[np.asarray(z["decomposition"][component])/z["decomposition"]["total"][0] for z in refs[k]]
            lineband(ax,np.arange(51)*.04,[v[:51] for v in y],palette[component],component.capitalize())
        ax.set(title=LAB[k],xlabel="Lag (s)",ylabel="Covariance / population variance",ylim=(None,1.38))
        ax.legend(frameon=False,fontsize=7,loc="upper right")
    scope(fig,whole+" · shaded bands span network seeds")
    save(fig,out,3,"covariance_decomposition")

    fig,axes=grid(2,len(conditions),(9.3,6.0))
    for i,k in enumerate(conditions):
        a,b=axes[i],axes[i+len(conditions)]
        current=[r for r in fit["paired"]["within_episode_shift"] if r["condition"]==k]
        for j,r in enumerate(current):
            shift=np.asarray(r["replicates"])
            a.plot([j-.13,j+.13],[r["original"],r["shifted_mean"]],color="#62747F",lw=1.)
            a.scatter(j-.13,r["original"],s=34,color=CT[k],zorder=4)
            a.scatter(np.full(len(shift),j+.13),shift,s=20,color="#798C94",alpha=.58,zorder=3)
        a.set(title=LAB[k],ylabel="d₂ (bits/s)",xticks=range(len(current)),
              xticklabels=[str(r["connectivity_seed"]) for r in current],xlabel="Paired delay seed")
        a.legend(handles=[Line2D([],[],marker="o",ls="",color=CT[k],label="Original"),
                          Line2D([],[],marker="o",ls="",color="#798C94",label="Episode shift")],
                 frameon=False,ncol=2,fontsize=7,loc="upper left")
        rows=refs[k]
        lineband(b,np.arange(51)*.04,[r["population"]["acf"][:51] for r in rows],CT[k],"Original")
        lineband(b,np.arange(51)*.04,[np.mean([s["population"]["acf"][:51] for s in r["controls"]
                 if s["kind"]=="shift"],axis=0) for r in rows],"#627F89","Episode shift","--")
        b.set(xlabel="Lag (s)",ylabel="Population normalized ACF",ylim=(None,1.28))
        b.legend(frameon=False,fontsize=7)
    scope(fig,whole+" · three independent shifts per network, fixed episode boundaries")
    save(fig,out,4,"circular_shift")

    if pulses:
        fig,axes=grid(2,2,(9.3,6.0))
        for i,key in enumerate(conditions):
            ax=axes[i]
            for state,color in (("up","#B56B19"),("down","#366C8B")):
                rows=[x for x in pulses if x["condition"]==key and x["pre_state"]==state]
                lineband(ax,np.arange(110)*.01-.1,[r["groups"]["untargeted_E"]["difference_hz"] for r in rows],
                         color,f"Pre-pulse {state.upper()} (n={len(rows)})")
            ax.axvline(0,c="#929EA3",lw=.8)
            ax.set(title=LAB[key],xlabel="Time from pulse (s)",ylabel="Untargeted E response (Hz/cell)")
            ax.legend(frameon=False,fontsize=7,loc="best")
        by={(r["context"],r["connectivity_seed"]):r["fit"]["d2"] for r in fs.values()
            if r["category"]=="original" and r["fit"]["d2"] is not None}
        for ax,k,l in ((axes[2],"integrated_250ms_spikes_per_cell","Extra spikes / untargeted E cell (250 ms)"),
                       (axes[3],"recovery_seconds","Untargeted E recovery time (s)")):
            for p in pulses:
                x=by[p["condition"],p["connectivity_seed"]]
                y=p["groups"]["untargeted_E"][k]
                m="o" if p["pre_state"]=="up" else "s"
                if y is not None:
                    ax.scatter(x,y,c=CT[p["condition"]],s=27,marker=m,alpha=.72)
                elif k=="recovery_seconds":
                    ax.scatter(x,1.,c=CT[p["condition"]],s=27,marker="^")
            ax.set(xlabel="Unperturbed d₂ (bits/s)",ylabel=l)
        scope(fig,whole+" · matched 10 ms pulse to 401 E cells, ranges across seeds")
        save(fig,out,6,"paired_pulse_response")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("results",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    r=a.results
    pulses=read(r/"pulses.json") if (r/"pulses.json").is_file() else None
    plot(read(r/"analysis.json"),read(r/"physiology.json"),read(r/"baseline_fits.json"),a.out,pulses)


if __name__=="__main__":
    main()

"""Render animal-level direction, state estimates and AR diagnostics."""
import argparse
from collections import defaultdict
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from plot_sleep_dandi import PALETTE, save

OLD = {"PNAS mouse V1":"Mouse", "000041":"Rat", "000552":"Mouse", "000061":"Rat", "hc-23":"Rat"}


def label(dataset, region):
    prefixes={"000056":"000056", "000939":"000939", "zenodo:19462601":"WT"}
    return f"{prefixes[dataset]} · {region}"


def directions(data,previous,out):
    new_rows=[x for x in data["animal_effects"] if x["context"]=="REST" and x["variant"]=="reference"]
    ids=sorted({(x["dataset"],x["region"]) for x in new_rows})
    subset=[x for x in previous if x["contrast"] in ("NREM_minus_Wake","REM_minus_NREM") and x["dataset"] in OLD]
    old_ids=list(dict.fromkeys((x["dataset"],x["region"],x["context"]) for x in subset))
    fig,axes=plt.subplots(1,2,figsize=(8.6,max(4.0,.24*(len(ids)+len(old_ids))+1.2)))
    fig.subplots_adjust(left=.18,right=.98,bottom=.10,top=.88,wspace=.60)
    labels=[label(*p) for p in ids]+[
        f"{OLD[d]} {r}"+(f" {c}" if c in ("PRE","POST") else "") for d,r,c in old_ids]
    all_effect=[x["effect"] for x in new_rows]+[v for x in subset for v in x["effect_by_replicate"].values()]
    mx=max(abs(x) for x in all_effect)*1.1
    for i,(ax,(title,short,long)) in enumerate(zip(axes,(("NREM − Wake","NW","NREM_minus_Wake"),
                                                        ("REM − NREM","RN","REM_minus_NREM")),strict=True)):
        groups=[[v["effect"] for v in new_rows if (v["dataset"],v["region"])==k and v["contrast"]==short] for k in ids]
        for d,r,c in old_ids:
            old=next((x for x in subset if (x["dataset"],x["region"],x["context"])==(d,r,c) and x["contrast"]==long),None)
            groups.append(list(old["effect_by_replicate"].values()) if old else [])
        ax.axvline(0,color=".55",lw=.7)
        ax.axhline(len(ids)-.5,color=".82",lw=.6)
        for y,points in enumerate(groups):
            if not points:
                continue
            new=y<len(ids)
            hollow=not new and old_ids[y-len(ids)][0]=="PNAS mouse V1"
            color="#135d77" if new else ".46"
            dy=np.linspace(-.18,.18,len(points)) if len(points)>1 else [0.]
            ax.scatter(points,y+np.asarray(dy),s=15,linewidth=.6,alpha=.7,
                       edgecolors=color,facecolors="white" if hollow else color,zorder=2)
            ax.scatter(np.mean(points),y,marker="|",s=95,color=color,linewidths=1.2,zorder=3)
        ax.set(yticks=range(len(labels)),yticklabels=labels,ylim=(len(labels)-.5,-.5),xlim=(-mx,mx),
               xlabel="Δd₂ (bits/s)",title=title)
        ax.text(0,1.075,chr(97+i),transform=ax.transAxes,fontsize=9,fontweight="bold")
    save(fig,out/"cross_region_directions")


def states(data,out):
    entries=[x for x in data["animal_state_d2"] if x["context"]=="REST"]
    groups=sorted({(x["dataset"],x["region"]) for x in entries})
    for ds in sorted({d for d,_ in groups}):
        regions=[r for d,r in groups if d==ds]
        for region in regions:
            values=defaultdict(dict)
            for x in entries:
                if x["dataset"]==ds and x["region"]==region:
                    values[x["subject"]][x["state"]]=x["d2"]
            order=sorted(values)
            if not order:
                continue
            fig,ax=plt.subplots(figsize=(max(4.9,.25*len(order)+1.1),2.65))
            fig.subplots_adjust(left=.13,right=.98,bottom=.28,top=.79)
            fig.text(.13,.925,label(ds,region),ha="left",va="center",fontsize=8.5)
            for j,id_ in enumerate(order):
                state=values[id_]
                placements=[(j+(i-1)*.19,s,state[s]) for i,s in enumerate(("Wake","NREM","REM")) if s in state]
                ax.plot([x for x,_,_ in placements],[y for _,_,y in placements],color=".73",lw=.6)
                for x,s,y in placements:
                    ax.scatter(x,y,c=PALETTE[s],s=20,zorder=3)
            ax.set(xticks=range(len(order)),xticklabels=order,ylabel="d₂ (bits/s)",xlim=(-.5,len(order)-.5))
            plt.setp(ax.get_xticklabels(),rotation=65,ha="right",rotation_mode="anchor",fontsize=7)
            fig.legend(handles=[Line2D([],[],marker="o",ls="",color=c,label=s) for s,c in PALETTE.items()],
                       loc="center right",bbox_to_anchor=(.98,.925),ncol=3,fontsize=7)
            key="WT" if ds.startswith("zenodo") else ds
            save(fig,out/f"state_d2_{key}_{region}")


def diagnostic(data,out):
    info=[x for x in data["diagnostics"] if x["context"]=="REST"]
    by_state=defaultdict(list)
    for x in info:
        by_state[x["state"]].append(x)
    fig,axs=plt.subplots(1,2,figsize=(7.9,2.9))
    fig.subplots_adjust(left=.105,right=.98,bottom=.20,top=.85,wspace=.36)
    axes=axs
    k=np.arange(1,21)
    for i,(state,color) in enumerate(PALETTE.items()):
        group=by_state[state]
        if not group:
            continue
        observed=np.median([x["observed_acf"] for x in group],axis=0)
        model=np.median([x["model_acf"] for x in group],axis=0)
        axes[0].plot(k,observed,color=color,lw=1)
        axes[0].plot(k,model,color=color,lw=.75,linestyle="--")
        tails=[x["residual_excess_kurtosis"] for x in group]
        axes[1].scatter(np.full(len(tails),i)+np.linspace(-.17,.17,len(tails)),tails,color=color,
                        s=8,alpha=.3)
        axes[1].scatter(i,np.median(tails),marker="_",color=color,s=80,lw=1.1)
    axes[0].set(xlabel="Lag (40 ms)",ylabel="Median autocorrelation",xticks=[1,5,10,15,20])
    axes[0].legend(handles=[Line2D([],[],color=".3",ls=mode,label=kind) for mode,kind in
                              (("-","Observed"),("--","AR(5)"))],fontsize=7)
    axes[1].set(xticks=range(3),xticklabels=list(PALETTE),ylabel="Residual excess kurtosis",yscale="symlog")
    for i,ax in enumerate(axes):
        ax.text(0,1.09,chr(97+i),transform=ax.transAxes,fontweight="bold",fontsize=9)
    save(fig,out/"ar_model_diagnostics")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("summary",type=Path)
    p.add_argument("historical",type=Path)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    x=json.loads(a.summary.read_text())
    prev=json.loads(a.historical.read_text())
    directions(x,prev,a.output)
    states(x,a.output)
    diagnostic(x,a.output)


if __name__=="__main__":
    main()

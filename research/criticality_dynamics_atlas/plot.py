"""Build the empirical atlas PDFs from the published lightweight summaries."""
from collections import defaultdict
from pathlib import Path
import argparse
import gzip
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np

from summarize import STATES, group_id

ROOT=Path(__file__).resolve().parent
COL=dict(Wake="#007484",NREM="#AC430C",REM="#71549F")
LAYERS={"historical_population":"historical","unit_matched":"180 s","NWB_matched":"three-state","PAG_1ms":"1 ms source"}
plt.rcParams.update({"font.size":8,"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
                    "pdf.fonttype":42,"axes.spines.top":False,"axes.spines.right":False,
                    "axes.labelcolor":"#172737","text.color":"#172737","axes.titleweight":"bold"})


def read(root,name):
    p=root/name
    return json.loads(gzip.decompress(p.read_bytes()) if p.suffix==".gz" else p.read_bytes())


def tag(g):
    a,b,c,d=g
    ctx={"natural":"", "baseline":"", "full_sleep":"", "general":"", "all":""}.get(d,d)
    return " ".join(x for x in (c,b,ctx,LAYERS.get(a,a)) if x)


def values(rows,name):
    return [r[name] for r in rows if r.get(name) is not None]


def animal_rows(rows,metric):
    g=defaultdict(list)
    for row in rows:
        v=row.get(metric)
        if v is not None:
            g[row["meta"]["subject"]].append(v)
    out=[]
    for x in g.values():
        x=np.asarray(x,dtype=float)
        n=np.isfinite(x).sum(axis=0)
        out.append(np.divide(np.nansum(x,axis=0),n,where=n>0,
                             out=np.full(np.shape(n),np.nan)))
    return out


def band(ax,x,curves,color,label):
    v=np.array(curves,dtype=float)
    if not len(v):
        return
    use=np.isfinite(v).any(axis=0)
    x=np.asarray(x)[use]
    v=v[:,use]
    ax.plot(x,np.nanmedian(v,axis=0),color=color,label=f"{label} (n={len(v)})",lw=1.4)
    if len(v)>2:
        a,b=np.nanquantile(v,[.25,.75],axis=0)
        ax.fill_between(x,a,b,color=color,alpha=.13,lw=0)


def layout(rows,cols,size):
    f,ax=plt.subplots(rows,cols,figsize=size,layout="constrained",squeeze=False)
    for i,a in enumerate(ax.flat):
        a.text(-.09,1.025,chr(65+i),transform=a.transAxes,weight="bold",size=10)
        a.set_axisbelow(True)
    return f,list(ax.flat)


def align_grid(f):
    if len(f.axes)<3:
        return
    cols=f.axes[0].get_subplotspec().get_gridspec().ncols
    if cols<3:
        return
    f.canvas.draw()
    pos=[ax.get_position() for ax in f.axes[:cols]]
    gutter=max(b.x0-a.x1 for a,b in zip(pos,pos[1:]))
    left=pos[0].x0
    width=(pos[-1].x1-left-(cols-1)*gutter)/cols
    f.set_layout_engine("none")
    for i,ax in enumerate(f.axes):
        y=ax.get_position()
        ax.set_position([left+(i%cols)*(width+gutter),y.y0,width,y.height])


def save(f,out,num,name):
    align_grid(f)
    f.savefig(out/f"{num:02d}_{name}.pdf")
    plt.close(f)


def points(ax,rows,x,y,label=None,marker="o",alpha=.62):
    for state in STATES:
        q=[r for r in rows if r["state"]==state and r.get(x) is not None and r.get(y) is not None]
        ax.scatter([r[x] for r in q],[r[y] for r in q],s=12,c=COL[state],alpha=alpha,
                   marker=marker,linewidths=0,label=state if label else None)
    ax.set(xlabel=x,ylabel=y)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results",type=Path,default=ROOT/"results")
    p.add_argument("--out",type=Path,default=ROOT/"figures")
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    r=a.results
    obs=read(r,"observations.json.gz")
    sens=read(r,"surrogates.json.gz")
    states=read(r,"animal_states.json")
    effects=read(r,"paired.json")
    assoc=read(r,"associations.json")
    null=read(r,"null_models.json")
    sim=read(r,"simultaneous.json")
    hypotheses=read(r,"hypothesis_evidence.json")
    gr=defaultdict(list)
    gs=defaultdict(list)
    for z in obs:
        gr[group_id(z)].append(z)
    for z in states:
        gs[tuple(z["group"])].append(z)
    cohort=[g for g in sorted(gs, key=lambda g:(g[2],g[1],g[3],g[0]))
            if g[0] in {"unit_matched","NWB_matched","PAG_1ms"} and len({v["animal"] for v in gs[g]})>=3]
    eq=[g for g in cohort if g[0]=="unit_matched"]
    example=next(g for g in eq if g[1]=="000939" and g[2]=="PoSub")
    adn=next(g for g in eq if g[2]=="ADN")
    frontal=next((g for g in eq if g[1]=="000041"),example)
    anchor=max((g for g in cohort if g[0]=="NWB_matched"),key=lambda g:len(gs[g]))
    by_state=lambda g,s:[x for x in gr[g] if x["meta"]["state"]==s]

    f,axs=layout(1,5,(15,6.6))
    for j,(ax,(metric,title)) in enumerate(zip(axs,(("d2","d₂ (bits/s)"),("tau_population","Population τ (s)"),
                      ("m_MR","MR m"),("pr_corr_n","Correlation PR/N"),("fano","Population Fano")))):
        for i,g in enumerate(cohort):
            for k,s in enumerate(STATES):
                v=values([x for x in gs[g] if x["state"]==s],metric)
                if v:
                    y=len(cohort)-1-i+(k-1)*.23
                    ax.scatter(v,np.full(len(v),y),s=10,color=COL[s],alpha=.52)
                    ax.plot(np.median(v),y,marker="|",ms=8,mew=1.8,color=COL[s])
                    ax.text(1.02+k*.062,len(cohort)-1-i,str(len(v)),transform=ax.get_yaxis_transform(),
                            size=6,ha="left",va="center",color=COL[s])
        ax.set(xlabel=title,ylim=(-.6,len(cohort)-.35),yticks=np.arange(len(cohort)),
               yticklabels=list(map(tag,cohort[::-1])) if j==0 else [])
        if metric in {"d2","tau_population","fano"}:
            ax.set_xscale("log")
    axs[-1].legend(handles=[Line2D([],[],marker="o",color=COL[s],ls="",label=s) for s in STATES],
                   loc="upper right",bbox_to_anchor=(1.08,1.12),frameon=False,ncol=3)
    save(f,a.out,1,"empirical_overview")

    f,axs=layout(2,2,(10.9,8))
    for ax,g in zip(axs[:2],(anchor,frontal)):
        for k,state in enumerate(STATES):
            data=[v[2] for x in by_state(g,state) for v in x.get("neurons",[]) if v[2] is not None]
            if len(data)>=2 and np.ptp(data)>0:
                w=ax.violinplot(np.log10(data),positions=[k],widths=.66,showmedians=True)
                for b in w["bodies"]:
                    b.set_facecolor(COL[state]); b.set_alpha(.38)
                for j in ("cbars","cmins","cmaxes","cmedians"):
                    w[j].set_color(COL[state])
            ax.scatter(np.full(len(values([q for q in gs[g] if q["state"]==state],"tau_population")),k+.22),
                       np.log10(values([q for q in gs[g] if q["state"]==state],"tau_population")),
                       s=14,color=COL[state],alpha=.55,marker="D")
        ax.set(title=tag(g),ylabel="log₁₀ τ (s)",xticks=range(3),xticklabels=STATES,xlim=(-.55,2.56))
    for ax,metric,title in zip(axs[2:],("tau_population","neuron_tau_same_units"),("Population τ","Same-cell median τ")):
        q=[z for z in effects if z["contrast"] in {"NREM_Wake","REM_NREM"} and z["metric"]==metric and tuple(z["group"]) in cohort]
        chosen=[g for g in cohort if any(tuple(z["group"])==g for z in q)]
        for z in q:
            yp=chosen.index(tuple(z["group"])) + (.17 if z["contrast"]=="REM_NREM" else -.17)
            color=COL[z["contrast"].split("_")[0]]
            ax.scatter(list(z["animals"].values()),np.full(z["n"],yp),s=9,color=color,alpha=.45)
            if z["ci"]:
                ax.plot(z["ci"],[yp,yp],color=color,lw=1.6)
            ax.scatter(z["mean"],yp,marker="D",s=16,c=color,zorder=3)
        ax.axvline(0,color="#77818C",lw=.9)
        ax.set(xlabel=f"Paired change in {title.lower()} (s)",title=title,yticks=range(len(chosen)),
               yticklabels=[tag(g) for g in chosen])
    axs[-1].legend(handles=[Line2D([],[],marker="D",color=COL[k],label=lab) for k,lab in
                            (("NREM","NREM − Wake"),("REM","REM − NREM"))],
                    loc="upper right",frameon=False,fontsize=7)
    save(f,a.out,2,"timescales")

    f,axs=layout(2,2,(10.6,8))
    for ax,g in zip(axs[:2],(example,frontal)):
        points(ax,gs[g],"tau_population","d2",True)
        ax.set(xscale="log",yscale="log",title=tag(g),xlabel="Population τ (s)",ylabel="d₂ (bits/s)")
        q=[z for z in null if z["order"]==1 and z["d2"] is not None]
        ax.plot([np.median([z["truth"]["tau"] for z in q if z["scenario"]==name]) for name in sorted({z["scenario"] for z in q})],
                [np.median([z["d2"] for z in q if z["scenario"]==name]) for name in sorted({z["scenario"] for z in q})],
                "--",color="#424E58",lw=1,label="Official AR(1) control")
        ax.legend(frameon=False,fontsize=7)
    points(axs[2],states,"tau_AR_longest","d2")
    axs[2].set(xscale="log",yscale="log",xlabel="Longest AR mode (s)",ylabel="d₂ (bits/s)",title="Estimator-shared across cohorts")
    all_del=defaultdict(dict)
    for q in effects:
        g=tuple(q["group"])
        if g in eq and q["contrast"]=="NREM_Wake" and q["metric"] in {"d2","tau_population"}:
            all_del[g][q["metric"]]=q["animals"]
    for g,m in all_del.items():
        x,y=m.get("tau_population",{}),m.get("d2",{})
        keys=sorted(x.keys()&y.keys())
        if keys:
            axs[3].scatter([x[k] for k in keys],[y[k] for k in keys],s=13,alpha=.65,label=tag(g))
    axs[3].axhline(0,color="#77818C",lw=.8); axs[3].axvline(0,color="#77818C",lw=.8)
    axs[3].set(xlabel="NREM − Wake population τ (s)",ylabel="NREM − Wake d₂ (bits/s)")
    axs[3].legend(frameon=False,fontsize=6,loc="upper left",bbox_to_anchor=(1,1))
    save(f,a.out,3,"distance_timescale")

    f,axs=layout(2,3,(13,7.5))
    for s in STATES:
        d=values([z for z in gs[anchor] if z["state"]==s],"m_MR")
        if d:
            axs[0].hist(d,bins=np.linspace(0,1,22),histtype="step",density=True,color=COL[s],label=f"{s} (n={len(d)})")
    axs[0].set(xlabel="MR m",ylabel="Animal-level density",title=tag(anchor)); axs[0].legend(frameon=False,fontsize=7)
    points(axs[1],gs[anchor],"m_MR","tau_population")
    axs[1].set(xlabel="MR m")
    points(axs[2],gs[anchor],"m_MR","d2")
    axs[2].set(xlabel="MR m")
    axs[1].set(ylabel="Population τ (s)",yscale="log",title=tag(anchor))
    axs[2].set(ylabel="d₂ (bits/s)",yscale="log",title=tag(anchor))
    for state in STATES:
        q=[z for z in by_state(anchor,state) if z.get("m_MR") is not None]
        if q:
            z=sorted(q,key=lambda z:z["m_MR"])[len(q)//2]
            k=np.arange(len(z["mr_coeff"]))
            axs[3].plot(k[1:]*.04,z["mr_coeff"][1:],color=COL[state],label=state)
            fit=z["mr_fit"]
            axs[3].plot(k[1:]*.04,fit["amp"]*fit["m"]**k[1:],"--",color=COL[state])
    axs[3].set(xlabel="Lag (s)",ylabel="Multistep slope",title="Observed slope and fitted exponential")
    axs[3].legend(frameon=False,fontsize=7)
    for j,fraction in enumerate((1.,.3,.1)):
        q=[z for z in null if z["scenario"].startswith("branch_") and z["truth"]["sampling"]==fraction and z["m_MR"] is not None]
        axs[4].scatter([z["truth"]["m"] for z in q],[z["m_MR"] for z in q],s=13,alpha=.66,label=f"{fraction:.0%}")
    axs[4].plot([.15,1],[.15,1],":",color="#555")
    axs[4].set(xlabel="Known m",ylabel="Recovered MR m"); axs[4].legend(frameon=False,fontsize=7)
    z=[i for i in gr[example] if i.get("m_MR") is not None and i.get("m_half") is not None and i.get("m_short") is not None]
    for y,label in (("m_short","0.4 s lags"),("m_half","Half units")):
        axs[5].scatter([i["m_MR"] for i in z],[i[y] for i in z],s=11,alpha=.7,label=label)
    axs[5].plot([.1,1],[.1,1],":",color="#555")
    axs[5].set(xlabel="MR m (2 s lags, full sample)",ylabel="MR m (sensitivity)")
    axs[5].legend(frameon=False,fontsize=7); save(f,a.out,4,"branching")

    f,axs=layout(2,3,(13,7.7))
    nt=min(16,min(z["units"] for z in gr[example]))
    for state in STATES:
        sub=by_state(example,state)
        q=[z|{"top_modes":z["eigenspectrum"][:nt]} for z in sub if z.get("eigenspectrum")]
        band(axs[0],np.arange(1,nt+1),animal_rows(q,"top_modes"),COL[state],state)
        q=[z|{"pr_curve":[z.get("fixed_n",{}).get(str(n),[None,None])[1] for n in (5,10,20,40)]}
           for z in sub]
        band(axs[1],[5,10,20,40],animal_rows(q,"pr_curve"),COL[state],state)
    axs[0].set(title=tag(example),xlabel="Eigenvalue rank",ylabel="Fraction of variance",yscale="log")
    axs[1].set(title=tag(example),xlabel="Fixed N",ylabel="Median raw PR/N")
    axs[0].set_ylim(top=axs[0].get_ylim()[1]*2.8)
    axs[0].legend(frameon=False,fontsize=7,loc="upper right")
    q=[z for z in gr[example] if z.get("pr_first_window_n") is not None]
    axs[2].scatter([z["pr_raw_n"] for z in q],[z["pr_first_window_n"] for z in q],s=11,color="#007484")
    axs[2].plot([0,1],[0,1],":",color="#555")
    axs[2].set(xlabel="Raw PR/N at 180 s",ylabel="Raw PR/N at 60 s")
    points(axs[3],gs[example],"pr_raw_n","pr_corr_n")
    axs[3].set(xlabel="Raw PR/N",ylabel="Correlation PR/N")
    points(axs[4],gs[example],"pr_fixed10","d2")
    axs[4].set(xlabel="Fixed 10 raw PR/N",ylabel="d₂ (bits/s)",yscale="log")
    points(axs[5],gs[example],"pr_fixed10","tau_integrated",True)
    axs[5].set(xlabel="Fixed 10 raw PR/N",ylabel="Integrated τ (s)",yscale="log")
    axs[5].legend(frameon=False,fontsize=7); save(f,a.out,5,"dimension")

    f,axs=layout(2,2,(10.8,7.7))
    sh=[z for z in sens if z["kind"]=="shift" and z["meta"]["dandiset"]==example[1]
        and z["meta"]["region"]==example[2] and z["meta"]["context"]==example[3]]
    for ax,state in zip(axs[:2],("Wake","NREM")):
        for k,col in (("auto","#007484"),("cross","#AC430C"),("total","#392C73")):
            q=animal_rows(by_state(example,state),"cov_"+k)
            u=animal_rows([z for z in sh if z["meta"]["state"]==state],"short_cov_"+k)
            if q:
                ax.plot(np.arange(51)*.04,np.median(np.asarray(q)[:,:51],axis=0),color=col,label=k)
            if u:
                ax.plot(np.arange(51)*.04,np.median(u,axis=0),"--",color=col)
        ax.set(title=f"{state} · {tag(example)}",xlabel="Lag (s)",ylabel="Covariance (counts² / bin)")
        ax.legend(frameon=False,fontsize=7)
    a0=read(r,"surrogate_effects.json")
    for i,kind in enumerate(("shift","permute","iaaft")):
        for j,state in enumerate(STATES[:2]):
            z=[v for v in a0 if tuple(v["group"])==example and v["kind"]==kind and v["state"]==state]
            for ax,col in ((axs[2],"cross_fraction0"),(axs[3],"d2")):
                x=values(z,col)
                if x:
                    y=i+(j-.5)*.24
                    ax.scatter(x,np.full(len(x),y),s=11,color=COL[state],alpha=.55)
                    ax.scatter(np.median(x),y,marker="D",s=21,c=COL[state])
    for ax,title in zip(axs[2:],("Δ cross variance fraction","Δ d₂ (bits/s)")):
        ax.axvline(0,c="#77818C",lw=.8)
        ax.set(xlabel=title+" relative to observed",yticks=range(3),yticklabels=["Shift","Permute","IAAFT"])
    save(f,a.out,6,"collective_covariance")

    f,axs=layout(2,3,(13.1,7.7))
    for ax,k,y in zip(axs[:4],("fano","variance","silent_fraction","off_duration"),
                      ("Population Fano","Variance (counts² / bin)","Silent occupancy","Silence duration p95 (s)")):
        points(ax,gs[example],"d2",k)
        ax.set(xlabel="d₂ (bits/s)",ylabel=y,title=tag(example))
        ax.xaxis.set_major_locator(MaxNLocator(5))
    for state in STATES:
        blocks=by_state(anchor,state)
        per=[]
        for z in blocks:
            q=np.asarray(z["off"]["hist"])
            if q.sum():
                per.append(z|{"off_norm":(q/q.sum()).tolist()})
        band(axs[4],np.arange(len(per[0]["off_norm"])),animal_rows(per,"off_norm"),COL[state],state) if per else None
    axs[4].set(xlabel="Silence run duration (s, histogram edge)",ylabel="Event fraction",
               xticks=[2,5,8,11,12],xticklabels=["0.08","0.32","1.28","8","16+"])
    axs[4].set_ylim(top=axs[4].get_ylim()[1]*1.5)
    axs[4].legend(frameon=False,fontsize=7,loc="upper right")
    q=[z for z in sens if z["kind"]=="active_only" and z["d2"] is not None]
    m={k["id"]:k for k in obs if k.get("id")}
    for state in STATES[:2]:
        u=[z for z in q if z["meta"]["state"]==state and m[z["id"].removesuffix("_active_only")]["d2"] is not None]
        axs[5].scatter([m[z["id"].removesuffix("_active_only")]["d2"] for z in u],[z["d2"] for z in u],
                       s=12,c=COL[state],alpha=.64,label=state)
    span=axs[5].get_xlim()
    if q:
        axs[5].plot(span,span,":",color="#77818C")
    axs[5].set(xlabel="All matched intervals d₂ (bits/s)",ylabel="Active-only d₂ (bits/s)",title="Frontal activity-conditioned sensitivity")
    axs[5].legend(frameon=False,fontsize=7)
    save(f,a.out,7,"fluctuations")

    f,axs=layout(2,3,(13.1,7.7))
    for state in STATES:
        q=by_state(anchor,state)
        if q:
            band(axs[0],np.arange(201)*.04,animal_rows(q,"acf"),COL[state],state)
            zz=[z|{"norm_psd":(np.asarray(z["psd"])/sum(z["psd"])).tolist()} for z in q if z.get("psd")]
            if zz:
                band(axs[1],zz[0]["freq"][1:],[z[1:] for z in animal_rows(zz,"norm_psd")],COL[state],state)
            for i,z in enumerate(q):
                if z.get("knee_hz"):
                    axs[2].scatter(STATES.index(state)+(i-len(q)/2)*.025,z["knee_hz"],s=12,c=COL[state],alpha=.7)
            for j,width in enumerate((.02,.04,.08,.20)):
                zz=[z|{"tau":z.get("tau_integrated") if width==.04 else
                       next((v.get("tau_integrated") for v in z.get("bins",[]) if v["bin_s"]==width),None)}
                    for z in q]
                vv=animal_rows(zz,"tau")
                if vv:
                    axs[3].scatter(np.full(len(vv),j+(STATES.index(state)-1)*.18),vv,s=10,c=COL[state],alpha=.44)
                    axs[3].plot(j+(STATES.index(state)-1)*.18,np.median(vv),marker="D",ms=5,c=COL[state])
    axs[0].set(xlabel="Lag (s)",ylabel="ACF",title=tag(anchor),ylim=(None,1.5))
    axs[0].legend(frameon=False,fontsize=7,loc="upper right")
    axs[1].set(xlabel="Frequency (Hz)",ylabel="Relative PSD",xscale="log",yscale="log")
    axs[2].set(xticks=range(3),xticklabels=STATES,ylabel="Qualified PSD knee (Hz)")
    axs[3].set(xticks=range(4),xticklabels=["20","40","80","200"],xlabel="Native bin (ms)",ylabel="Integrated τ (s)")
    points(axs[4],gs[example],"d2","d2_AR10",True)
    axs[4].set(xlabel="40 ms AR(5) d₂",ylabel="40 ms AR(10) d₂",xscale="log",yscale="log")
    q=[z for z in sens if z["kind"]=="multiscale" and z["meta"]["dandiset"]==example[1]
       and z["meta"]["region"]==example[2] and z["meta"]["context"]==example[3]]
    for i,width in enumerate((.04,.08,.20)):
        items0=gr[example] if width==.04 else [z for z in q if z["bin_seconds"]==width]
        for j,s in enumerate(STATES[:2]):
            yy=values([z for z in items0 if z["meta"]["state"]==s],"d2")
            if yy:
                axs[5].scatter(np.full(len(yy),i+(j-.5)*.22),yy,s=12,c=COL[s],alpha=.54)
    axs[5].set(xticks=range(3),xticklabels=["40","80","200"],xlabel="Archive-derived bin (ms)",ylabel="AR(5) d₂ (bits/s)",yscale="log")
    save(f,a.out,8,"multiscale")

    f,axs=layout(1,3,(13.5,6.5))
    fk=("tau_integrated","pr_corr_n","cross_fraction0","fano","d2")
    for ax,contrast in zip(axs[:2],("NREM_Wake","REM_NREM")):
        q={(tuple(z["group"]),z["metric"]):z for z in effects if z["contrast"]==contrast}
        for i,g in enumerate(cohort):
            for j,k in enumerate(fk):
                z=q.get((g,k))
                if z:
                    signed=(2*z["positive"]-z["n"])/z["n"]
                    ax.scatter(j,i,c=[signed],cmap="RdBu_r",vmin=-1,vmax=1,marker="s",s=275)
                    ax.text(j,i,str(z["n"]),ha="center",va="center",fontsize=6,
                            color="white" if abs(signed)>.63 else "#172737")
        ax.set(xlim=(-.6,4.6),ylim=(-.6,len(cohort)-.4),xticks=range(5),
               xticklabels=["τ int.","PR/N","Cross","Fano","d₂"],yticks=range(len(cohort)),
               yticklabels=list(map(tag,cohort)) if ax==axs[0] else [],title=contrast.replace("_"," − "))
    zz=[x for x in sim if x.get("tau_integrated") is not None and
        (x["key"][1]=="000061" and x["key"][0]=="NWB_matched" or
         x["key"][1]=="zenodo:19462601" and x["key"][0]=="unit_matched")]
    eligible=defaultdict(set)
    for z in zz:
        k=z["key"]
        eligible[k[1],z["region_a"],z["region_b"],k[4]].add(k[2])
    pairs=sorted(q for q,n in eligible.items() if len(n)>=3)

    for y,pair in enumerate(pairs):
        for j,s in enumerate(STATES):
            gg=defaultdict(list)
            for z in zz:
                if (z["key"][1],z["region_a"],z["region_b"],z["key"][4])==pair and z["key"][5]==s:
                    gg[z["key"][2]].append(z["tau_integrated"])
            if gg:
                d=[np.mean(v) for v in gg.values()]
                sub=y+(j-1)*.21
                axs[2].scatter(d,np.full(len(d),sub),s=14,c=COL[s],alpha=.61)
                axs[2].scatter(np.median(d),sub,marker="D",s=21,c=COL[s])
    axs[2].axvline(0,c="#77818C",lw=.8)
    axs[2].set(xlabel="Same-time region A − B τ (s)",title="Exact simultaneous windows",
               yticks=range(len(pairs)),yticklabels=[f"{s} {v}−{w} {d}" for s,v,w,d in pairs])
    save(f,a.out,9,"region_state")

    f,axs=layout(2,3,(13.7,7.8))
    aa=[z for z in null if z["order"]==1 and z["d2"] is not None]
    for ax,k,l in zip(axs[:3],("tau","var","chi"),("Exact decay time (s)","Variance at unit innovation","Static response χ")):
        ax.scatter([z["truth"][k] for z in aa],[z["d2"] for z in aa],s=16,color="#007484",alpha=.65)
        ax.set(xlabel=l,ylabel="Official AR(1) d₂",xscale="log",yscale="log")
    for j,wgt in enumerate((.05,.25,.65)):
        z=[q for q in null if q["scenario"].startswith("mix_") and q["truth"]["weight"]==wgt and q["d2"] is not None]
        axs[3].scatter([q["truth"]["long_tau"] for q in z],[q["d2"] for q in z],s=15,alpha=.56,label=f"Slow weight {wgt:g}")
    axs[3].set(xscale="log",yscale="log",xlabel="Identical constituent τ within each column (s)",ylabel="AR(5) d₂")
    axs[3].legend(frameon=False,fontsize=7)
    z=[q for q in null if q["scenario"].startswith("branch_") and q["d2"] is not None]
    for h in (1.,.3,.1):
        zz=[q for q in z if q["truth"]["sampling"]==h]
        axs[4].scatter([q["truth"]["m"] for q in zz],[q["d2"] for q in zz],s=16,alpha=.62,label=f"{h:.0%}")
    axs[4].set(xlabel="True branching m",ylabel="AR(5) d₂",yscale="log"); axs[4].legend(frameon=False,fontsize=7)
    gates={"independent_poisson":"Poisson", "common_poisson":"Slow\ndrive","latent_ungated":"Ungated",
           "gate_0.15_3":"15%\n0.12 s","gate_0.15_15":"15%\n0.60 s",
           "gate_0.35_3":"35%\n0.12 s","gate_0.35_15":"35%\n0.60 s"}
    for i,g in enumerate(gates):
        q=[z for z in null if z["scenario"]==g and z["d2"] is not None]
        axs[5].scatter(np.full(len(q),i),values(q,"d2"),s=14,color="#AC430C",alpha=.6)
    axs[5].set(xticks=range(len(gates)),xticklabels=gates.values(),yscale="log",ylabel="AR(5) d₂",
               xlabel="Control or gate (OFF fraction, mean duration)")
    save(f,a.out,10,"theoretical_nulls")

    f,axs=layout(1,2,(12.7,7.5))
    for ax,g in zip(axs,(example,adn)):
        groupsel={z["y"]:z for z in assoc if tuple(z["group"])==g and z["x"]=="d2" and z["level"]=="delta:NREM_Wake"}
        shown=["tau_AR_longest","tau_population","tau_integrated","neuron_tau_median","m_MR","pr_corr_n",
               "cross_fraction0","fano","silent_fraction","off_duration","knee_hz","persist1s"]
        for i,k in enumerate(shown):
            z=groupsel.get(k)
            if z:
                color={"A":"#637078","B":"#007484","C":"#AC430C"}[z["type"]]
                if z["ci"]:
                    ax.plot(z["ci"],[i,i],c=color,lw=1.7)
                ax.scatter(z["rho"],i,s=20,c=color,zorder=3)
                ax.text(1.09,i,f"n={z['n_animals']}",ha="left",va="center",fontsize=7)
        ax.axvline(0,c="#77818C",lw=.8)
        ax.set(yticks=range(len(shown)),yticklabels=[k.replace("_"," ") for k in shown],
               ylim=(len(shown)-.5,-.5),xlim=(-1.05,1.34),xlabel="Spearman ρ for paired state changes",title=tag(g))
    save(f,a.out,11,"association_matrix")

    f,ax=plt.subplots(figsize=(11,4.1),layout="constrained")
    ax.set_axis_off()
    x=(.02,.09,.33)
    for i,h in enumerate(("Test","Judgment","Observed result")):
        ax.text(x[i],.92,h,transform=ax.transAxes,weight="bold",fontsize=11)
    for i,row in enumerate(hypotheses):
        y=.81-i*.111
        ax.axhline(y+.052,c="#CDD6D9",lw=.6)
        for j,k in enumerate(("hypothesis","judgment","evidence")):
            ax.text(x[j],y,row[k],transform=ax.transAxes,fontsize=9,va="center")
    save(f,a.out,12,"hypothesis_evaluation")


if __name__=="__main__":
    main()

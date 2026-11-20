"""Draw sleep results from verified MATLAB fits and animal-level summaries."""

from collections import defaultdict
from pathlib import Path
import argparse
import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from summarize_sleep_dandi import event_summary, state_summary, valid


mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans","sans-serif"],
                     "font.size":8,"svg.fonttype":"none","pdf.fonttype":42,
                     "axes.spines.right":False,"axes.spines.top":False,"axes.linewidth":0.7,"legend.frameon":False})
PALETTE = {"Wake":"#91604a","NREM":"#135d77","REM":"#705495"}
STATES = ("Wake","NREM","REM")


def save(fig, base):
    base.parent.mkdir(parents=True,exist_ok=True)
    if len(fig.axes)>1 and os.getenv("NATURE_FIGURE_AUDIT"):
        from audit_panel_alignment import require_matplotlib_panel_alignment
        out = Path(os.environ.get("NATURE_FIGURE_AUDIT_DIR",str(base.parent)))
        out.mkdir(parents=True,exist_ok=True)
        require_matplotlib_panel_alignment(fig,json_out=str(out/(base.name+".alignment.json")),
                                           tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True)
    for suffix in ("pdf","svg","png"):
        fig.savefig(base.with_suffix("."+suffix),dpi=600 if suffix=="png" else None,facecolor="white")
    plt.close(fig)


def pilot(data, out):
    rows = {r["record"]["state"]:r for r in data if r["record"]["variant"]=="reference" and valid(r["fit"])}
    assert {"Wake","NREM"}<=set(rows)
    fig,axes = plt.subplots(1,2,figsize=(7.15,2.65))
    fig.subplots_adjust(left=.105,right=.97,bottom=.20,top=.91,wspace=.42)
    use = [s for s in STATES if s in rows]
    for i,s in enumerate(use):
        row=rows[s]
        fit=row["fit"]
        axes[0].scatter(i,fit["d2"],color=PALETTE[s],s=35,zorder=3)
        if valid(row["long_reach"]):
            axes[0].scatter(i,row["long_reach"]["d2"],facecolors="white",edgecolors=PALETTE[s],s=33,marker="D")
        axes[1].plot(np.arange(1,21),fit["residual_acf"],color=PALETTE[s],linewidth=1.2,label=s)
    axes[0].set(xticks=range(len(use)),xticklabels=use,ylabel="d₂ (bits/s)",title="ACC  ·  BWRat17")
    axes[0].set_xlim(-.45,len(use)-.55)
    axes[1].axhline(0,color=".65",linewidth=.75)
    axes[1].set(xlabel="Lag (40 ms bins)",ylabel="Residual autocorrelation",xticks=[1,5,10,15,20])
    axes[1].legend(loc="best",fontsize=7)
    save(fig,out/"pilot_acc")


def cohort(rows,out):
    summary=state_summary(rows)
    paired=[r for r in summary["per_animal_dataset"] if r["variant"]=="reference" and r["context"]=="HOME"]
    if not paired:
        return
    subjects=sorted({r["subject"] for r in paired})
    nrem=[r for r in paired if r["contrast"]=="NREM_minus_Wake"]
    rem=[r for r in paired if r["contrast"]=="REM_minus_NREM"]
    fig,ax=plt.subplots(figsize=(7.15,3.12),layout="constrained")
    ax.axhline(0,color=".55",linewidth=.75)
    for diff,label,color,marker,offset in ((nrem,"NREM − Wake","#135d77","o",-.11),
                                          (rem,"REM − NREM","#705495","D",.11)):
        if diff:
            ax.scatter([subjects.index(r["subject"])+offset for r in diff],[r["effect"] for r in diff],
                       s=23,color=color,marker=marker,label=label,zorder=3)
    ax.set(xticks=range(len(subjects)),xticklabels=subjects,ylabel="Mean paired Δd₂ (bits/s)")
    ax.tick_params(axis="x",labelrotation=42,labelsize=7)
    ax.legend(loc="best")
    save(fig,out/"within_region")


def cross_region(rows,out):
    records=[r for r in state_summary(rows)["paired_cross_region_sessions"]
             if r["variant"]=="reference" and r["contrast"]=="NREM_minus_Wake"]
    if not records:
        return
    cells=defaultdict(list)
    for r in records: cells[r["subject"],r["context"]].append(r)
    fig,axes=plt.subplots(1,2,figsize=(7.15,2.85),sharey=True)
    fig.subplots_adjust(left=.105,right=.97,bottom=.18,top=.88,wspace=.27)
    colors=["#135d77","#86543d","#705495","#457653"]
    subjects=sorted({r["subject"] for r in records})
    for ax,context in zip(axes,("PRE","POST"),strict=True):
        ax.axhline(0,color=".55",linewidth=.75)
        for s,color in zip(subjects,colors):
            records=cells.get((s,context),[])
            if records:
                ax.plot([0,1],[np.mean([p[k] for p in records]) for k in ("ca1","bla")],
                        color=color,marker="o",markersize=4,label=s,lw=1.2)
        ax.set(xticks=[0,1],xticklabels=["CA1","BLA"],xlim=[-.35,1.35],title=context)
    axes[0].set_ylabel("Mean paired Δd₂ (bits/s)")
    handles,labels=axes[0].get_legend_handles_labels()
    if labels:
        fig.legend(handles,labels,loc="upper center",ncol=len(labels),fontsize=7,bbox_to_anchor=(.55,.99))
    save(fig,out/"paired_ca1_bla")


def swr_figure(rows,seed,out):
    details=event_summary(rows,seed)
    groups=defaultdict(list)
    for r in details:
        if r["status"]=="DESCRIPTIVE":
            groups[r["subject"]].append(r["matched_difference"]["median"])
    if not groups:
        return
    subjects=sorted(groups)
    fig,ax=plt.subplots(figsize=(7.15,2.9),layout="constrained")
    ax.axhline(0,color=".55",linewidth=.75,zorder=0)
    for i,animal in enumerate(subjects):
        values=groups[animal]
        ax.scatter([i]*len(values),values,s=12,color="#91b1bc",zorder=2)
        ax.scatter(i,np.mean(values),s=26,color="#135d77",zorder=3)
    ax.set(xticks=range(len(subjects)),xticklabels=subjects,ylabel="Matched SWR-rich − poor Δd₂ (bits/s)")
    ax.tick_params(axis="x",labelrotation=42,labelsize=7)
    save(fig,out/"swr_matched")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("kind",choices=["pilot","cohort","cross","events"])
    p.add_argument("input",type=Path)
    p.add_argument("output",type=Path)
    o=p.parse_args()
    fits=json.loads(o.input.read_text())
    rows=fits["records"]
    if dict is type(rows): rows=[rows]
    if o.kind=="events":
        swr_figure(rows,fits["settings"]["seed"],o.output)
    else:
        {"pilot":pilot,"cohort":cohort,"cross":cross_region}[o.kind](rows,o.output)

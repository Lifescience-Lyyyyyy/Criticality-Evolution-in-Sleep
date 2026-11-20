"""Render Phase 9 sleep evidence and the independently verified reference cohorts."""

from collections import defaultdict
from pathlib import Path
from statistics import mean
import argparse
import json

import matplotlib.pyplot as plt
import numpy as np

from plot_sleep_dandi import PALETTE, save
from summarize_sleep_dandi import valid

CONTRASTS = (("NREM_minus_Wake", "NREM − Wake"), ("REM_minus_NREM", "REM − NREM"))
CONTEXTS = ("PRE", "POST")
COLORS = {"PRE": "#135d77", "POST": "#9a5832"}


def paired_mec(rows, output):
    session = defaultdict(dict)
    for row in rows:
        r = row["record"]
        if r["variant"] == "reference" and r["dandiset"] == "hc-23" and valid(row["fit"]):
            session[r["subject"], r["session"]][r["state"]] = row["fit"]["d2"]
    animals = defaultdict(list)
    for (subject, _), values in session.items():
        if {"NREM", "REM"} <= set(values):
            animals[subject].append([values["NREM"], values["REM"]])
    if not animals:
        raise ValueError("No verified MEC sleep pairs")
    fig, ax = plt.subplots(figsize=(5.8,3.15),layout="constrained")
    subjects = sorted(animals)
    slots = np.arange(len(subjects))
    values = np.array([np.mean(animals[s],axis=0) for s in subjects])
    for x, y in zip(slots, values, strict=True):
        ax.plot(x + np.array([-.16,.16]),y,color=".64",linewidth=.8)
    for i, (state, marker) in enumerate((("NREM","s"),("REM","D"))):
        ax.scatter(slots + (-.16 if i==0 else .16),values[:,i],s=35,c=PALETTE[state],marker=marker,label=state,zorder=3)
    ax.set(xticks=slots,xticklabels=subjects,ylabel="Mean paired d₂ (bits/s)")
    ax.legend(loc="best")
    save(fig,output/"mec_paired_nrem_rem")


def ordered_labels(groups):
    labels = {("PNAS mouse V1","V1"):"Mouse V1",
              ("000041","Frontal cortex"):"Rat frontal cortex", ("000552","CA1"):"Mouse CA1",
              ("000061","CA1"):"Rat CA1", ("000061","BLA"):"Rat BLA",
              ("hc-23","MEC II/III"):"Rat MEC II/III"}
    rank = {key:i for i,key in enumerate(labels)}
    out = []
    for row in sorted(groups,key=lambda r:(rank[r["dataset"],r["region"]],r["context"]=="POST")):
        label = labels[row["dataset"],row["region"]]
        if row["context"] in CONTEXTS:
            label += " " + row["context"]
        out.append((row,label))
    return out


def unified(groups, output):
    fig, axes = plt.subplots(1,2,figsize=(8.4,4.15))
    fig.subplots_adjust(left=.19,right=.982,bottom=.16,top=.88,wspace=.82)
    panels = [[(r,l) for r,l in ordered_labels(groups) if r["contrast"] == contrast] for contrast,_ in CONTRASTS]
    maximum = max(len(p) for p in panels)
    vals = np.array([v for p in panels for row,_ in p for v in row["effect_by_replicate"].values()])
    span = max(abs(vals))*1.1
    for i,(ax,panel,(_,title)) in enumerate(zip(axes,panels,CONTRASTS,strict=True)):
        ax.axvline(0,color=".68",linewidth=.75)
        ax.set(xlim=(-span,span),ylim=(maximum-.4,-.6),
               yticks=range(len(panel)),yticklabels=[label for _,label in panel],xlabel="Δd₂ (bits/s)")
        ax.set_title(title,pad=8)
        ax.text(0,1.08,chr(97+i),transform=ax.transAxes,fontweight="bold",fontsize=8)
        for j,(row,_) in enumerate(panel):
            points = list(row["effect_by_replicate"].values())
            jitter = np.linspace(-.15,.15,len(points)) if len(points)>1 else [0]
            edge = "#135d77" if row["dataset"]=="hc-23" else ".32"
            tone = "white" if row["replicate"]=="session" else edge
            ax.scatter(points,j+np.array(jitter),facecolors=tone,edgecolors=edge,linewidth=.65,s=15,alpha=.7,zorder=2)
            ax.scatter(row["effect_mean"],j,marker="|",s=135,color=edge,lw=1.3,zorder=3)
    save(fig,output/"cross_region_sleep_effects")


def simultaneous(cross, output):
    fig, axes = plt.subplots(1,2,figsize=(7.2,2.9),sharey=True)
    fig.subplots_adjust(left=.115,right=.982,bottom=.21,top=.83,wspace=.27)
    for i,(ax,(contrast,name)) in enumerate(zip(axes,CONTRASTS,strict=True)):
        pairs = defaultdict(list)
        for row in cross["paired_cross_region_sessions"]:
            if row["variant"]=="reference" and row["contrast"]==contrast:
                pairs[row["subject"],row["context"]].append(row["ca1_minus_bla"])
        subjects = sorted({subject for subject,_ in pairs})
        slots = np.arange(len(subjects))
        ax.axhline(0,color=".67",linewidth=.75)
        for ctx,shift,marker in (("PRE",-.12,"o"),("POST",.12,"D")):
            x = [k+shift for k,s in zip(slots,subjects,strict=True) if (s,ctx) in pairs]
            y = [mean(pairs[s,ctx]) for s in subjects if (s,ctx) in pairs]
            ax.scatter(x,y,s=32,color=COLORS[ctx],marker=marker,label=ctx,zorder=3)
        ax.set(xticks=slots,xticklabels=subjects,title=name)
        ax.text(0,1.11,chr(97+i),transform=ax.transAxes,fontweight="bold",fontsize=8)
    axes[0].set_ylabel("CA1 − BLA state effect (bits/s)")
    axes[1].legend(loc="best",fontsize=7)
    save(fig,output/"phase8_same_time_ca1_bla")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mec",type=Path)
    p.add_argument("comparison",type=Path)
    p.add_argument("cross",type=Path)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    paired_mec(json.loads(a.mec.read_text())["records"],a.output)
    unified(json.loads(a.comparison.read_text()),a.output)
    simultaneous(json.loads(a.cross.read_text()),a.output)

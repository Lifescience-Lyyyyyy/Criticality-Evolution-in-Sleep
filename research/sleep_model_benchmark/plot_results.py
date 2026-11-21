"""Recreate single-panel validation, observation and multiscale figures."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np

DIR = Path(__file__).parent
COLOR = "#17676F"
ACCENT = "#8A4562"
GREY = "#757575"

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 9,
                     "axes.spines.top":False,"axes.spines.right":False,
                     "svg.fonttype":"none","svg.hashsalt":"sleep-model-benchmark","pdf.fonttype":42})


def export(fig,out,stem,preview):
    fig.savefig(out/f"{stem}.pdf",metadata={"CreationDate":None,"ModDate":None})
    svg=out/f"{stem}.svg"
    fig.savefig(svg,metadata={"Date":None})
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines())+"\n",encoding="utf-8")
    if preview:
        fig.savefig(preview/f"{stem}.png",dpi=300)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,default=DIR/"figures")
    p.add_argument("--preview",type=Path)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    if a.preview:
        a.preview.mkdir(parents=True,exist_ok=True)
    d=json.loads((DIR/"results/matlab_benchmark.json").read_text())
    original=json.loads((DIR/"results/original_validation.json").read_text())["results"]["rat_frontal"]
    fig,ax=plt.subplots(figsize=(5.25,2.65),layout="constrained")
    metrics=["up_mean_s","down_mean_s","up_cv","down_cv"]
    ax.axvspan(-2,2,color="#F2F4F4",label="Published mean ± 2 recording SD")
    ax.axvline(0,color=GREY,lw=.8)
    for i,key in enumerate(metrics):
        target,sd=original["reference"][key]
        vals=[(r[key]-target)/sd for r in original["validation"]]
        fill="white" if key==original["calibration_metric"] else COLOR
        ax.scatter(vals,np.full(3,i)+[-.18,0,.18],edgecolors=COLOR,facecolors=fill,s=27,zorder=4)
    ax.set(xlabel="Deviation from published mean (recording SD)",yticks=range(4),
           yticklabels=["UP mean","DOWN mean (time calibration)","UP CV","DOWN CV"],xlim=(-2.55,2.55))
    ax.invert_yaxis()
    ax.legend(loc="lower right",bbox_to_anchor=(1,1.025),frameon=False,fontsize=8)
    export(fig,a.out,"physiological_validation",a.preview)

    fig,ax=plt.subplots(figsize=(5.3,2.75),layout="constrained")
    labels={"original_poisson":"Poisson","original_shared_gamma":"Shared Gamma–Poisson",
            "original_half_rate":"Poisson at 0.5× rate","original_double_rate":"Poisson at 2× rate"}
    for i,(context,label) in enumerate(labels.items()):
        rows=sorted((r for r in d["rows"] if r["record"]["spec"]=="primary" and r["record"]["context"]==context),
                    key=lambda r:r["record"]["seed"])
        for j,r in enumerate(rows):
            y=i+(j-1)*.14
            ax.plot([f["d2"] for f in r["fits"]],[y,y],c="#BFC7C7",lw=.8,zorder=1)
            for f,color in zip(r["fits"],(COLOR,ACCENT),strict=True):
                ax.scatter(f["d2"],y,color=color,s=21,zorder=2,
                    label=f'AR({f["order"]})' if i==j==0 else None)
    ax.set(xlabel="d₂ (bits/s)",yticks=range(4),yticklabels=list(labels.values()),xlim=(1.35,3.6))
    ax.invert_yaxis()
    ax.legend(loc="upper right",frameon=False,fontsize=8)
    export(fig,a.out,"observation_sensitivity",a.preview)

    fig,ax=plt.subplots(figsize=(5.25,3.0),layout="constrained")
    model=[r for r in d["rows"] if r["record"]["context"]=="original_poisson"]
    policies=[("Fixed AR(5)",COLOR,{20:5,40:5,80:5,200:5}),
              ("History near 200 ms",ACCENT,{20:10,40:5,80:2,200:1})]
    for label,color,grid in policies:
        centers=[]
        for width,order in grid.items():
            values=[f["d2"] for r in model if round(1000*r["record"]["bin_seconds"])==width
                    and r["record"]["spec"]!="ms80_high" for i,f in enumerate(r["fits"])
                    if f["order"]==order and f["pass_screen"]]
            ax.scatter([width]*len(values),values,color=color,alpha=.7,s=18)
            centers.append(np.median(values))
        ax.plot(list(grid),centers,marker="o",lw=1.25,ms=3,color=color,label=label)
    alt=[r["fits"][0]["d2"] for r in model if r["record"]["spec"]=="ms80_high"]
    ax.scatter([80]*len(alt),alt,c=ACCENT,s=17,marker="x",label="AR(3) at 80 ms, 240 ms history")
    ax.set(xscale="log",xticks=[20,40,80,200],xticklabels=["20","40","80","200"],
           xlabel="Count bin (ms)",ylabel="d₂ (bits/s)",ylim=(0,6.35))
    ax.minorticks_off()
    ax.legend(loc="upper right",frameon=False,fontsize=8)
    export(fig,a.out,"multiscale_sensitivity",a.preview)


if __name__=="__main__":
    main()

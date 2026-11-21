"""Audit biological identities, repeat measures, provenance overlap and saved aggregates."""
from collections import defaultdict
from itertools import combinations
from statistics import mean
from hashlib import sha256
import math
import re

from build_atlas_data import ROOT, aggregate, build, read, write
from catalog import BASE, PAIRS, V1_FILE
from validate_atlas_data import verify_summaries
from build_pairwise_fixture import OUT
from select_preliminary_sample import BASELINE


def token(s):
    return re.sub(r"[^a-z0-9]","",str(s).lower())


def digests(v):
    return set(v.values()) if isinstance(v,dict) else {v}


def equal(a,b):
    return math.isclose(a,b,abs_tol=1e-12,rel_tol=0)


def run():
    regions,datasets,_,index=build()
    all_cohorts={c["id"]:c for choices in index.values() for c in choices}
    # This verifier compares against independent Phase 8–10 animal/session summaries.
    verify_summaries(all_cohorts, ROOT/"website/public/data")
    identities=defaultdict(lambda:defaultdict(set))
    assets=defaultdict(dict)
    summaries=[]
    for c in all_cohorts.values():
        species,region,dataset,ctx=c["id"].split("|")
        samples=c["records"]
        assert species==datasets[dataset]["species"]
        subjects={r["subject"] for r in samples}
        assert len(samples)==len({r["id"] for r in samples})
        rows={}
        for contrast,(a,b) in PAIRS.items():
            for model in ("reference","AR10"):
                pairs=defaultdict(list)
                for r in samples:
                    fits=r["fits"].get(model,{})
                    if a not in fits or b not in fits:
                        continue
                    subject=r["subject"] or r["session"]
                    pairs[subject].append(fits[b]["value"]-fits[a]["value"])
                saved=aggregate(samples,model,[a,b])
                counts={x["id"]:x for x in saved}
                assert set(pairs)==set(counts)
                for subject,vals in pairs.items():
                    assert counts[subject]["sessions"]==len(vals) and equal(counts[subject]["effect"],mean(vals))
                exported=c["effects"].get(model,{}).get(contrast)
                if exported:
                    assert exported["n"]==len(pairs) and exported["sessions"]==sum(map(len,pairs.values()))
                    assert not pairs or equal(exported["mean"],mean(mean(d) for d in pairs.values()))
            rows[contrast]=c["effects"]["reference"][contrast]
        assert c["animals"]==(None if None in subjects else len(subjects))
        for r in samples:
            subject=r["subject"]
            ids=assets[dataset]
            a=ids.setdefault(r["id"],dict(subject=subject,session=r["session"],asset=r["asset"],hashes=digests(r["digest"])))
            assert a["subject"]==subject and a["session"]==r["session"] and a["asset"]==r["asset"]
            if subject:
                identities[dataset][subject].add((region,ctx,r["id"]))
        summaries.append({"id":c["id"],"species":species,"dataset":dataset,"context":ctx,"region":region,
                          "animals":c["animals"],"sessions":len(samples),"states":c["states"],"effects":rows,
                          "count_unit":"session" if None in subjects else "animal"})
    v1=all_cohorts["mouse|V1|PNAS-V1|PUBLISHED"]
    assert v1["animals"] is None and v1["sessions"]==13 and v1["effects"]["reference"]["RN"]["n"]==0
    assert len(read(V1_FILE)["paired_nrem_minus_wake"])==13
    repeat=[]
    for ds, counts in identities.items():
        if ds=="000061":
            pre={s for s,k in counts.items() if any(c=="PRE" for _,c,_ in k)}
            post={s for s,k in counts.items() if any(c=="POST" for _,c,_ in k)}
            assert pre==post=={"Rat08","Rat10","Rat11"}
        if ds=="000056":
            shared=set.intersection(*({s for s, observations in counts.items()
                                        if any(r==region and c=="REST" for r,c,_ in observations)}
                                      for region in ("ADN","PoSub")))
            assert shared=={"Mouse24","Mouse25","Mouse28"}
        repeat.append(dict(dataset=ds,animals=len(counts),
                           multiple_sessions=sum(len({p for _,_,p in k})>1 for k in counts.values()),
                           multiple_regions=sum(len({r for r,_,_ in k if r!="Frontal"})>1 for k in counts.values())))
    publications={ds:dict(species=o["species"],title=o["title"],paper=o["paper"],
                          published_version=o["version"],release=o["url"],
                          subject_ids=sorted(identities[ds]),unique_assets=len(assets[ds]))
                  for ds,o in datasets.items()}
    checks=[]
    for a,b in combinations(sorted(assets),2):
        x,y=assets[a],assets[b]
        id_overlap=sorted(set(x)&set(y))
        fingerprint={h for r in x.values() for h in r["hashes"] if h}
        digest_overlap=sorted(fingerprint & {h for r in y.values() for h in r["hashes"] if h})
        sig=lambda c:{(token(r["subject"]),token(r["session"])) for r in c.values() if r["subject"]}
        checks.append(dict(datasets=[a,b],species_match=datasets[a]["species"]==datasets[b]["species"],
                           shared_publication_DOI=datasets[a]["paper"]==datasets[b]["paper"],
                           shared_asset_id=id_overlap,shared_digest=digest_overlap,
                           matching_subject_session=sorted(sig(x)&sig(y))))
        assert not (id_overlap or digest_overlap or sig(x)&sig(y)),(a,b)
    receipt=dict(baseline=BASELINE,source_metadata=publications,repeat_measures=repeat,inter_dataset_checks=checks,
                 subjects_ACROSS_datasets_aggregated=False,cohorts=summaries,
                 notes=[
                    "同场的状态值先配对，重复记录在动物内等权。比较独立发布来源的归档资产 ID、文件摘要、归一化动物与会话键、DOI 及记录协议。",
                    "000041、000056、000061、000552 有研究团队或公共归档来源关联。动物/会话编码、记录部位或行为协议不同，未找到同一原始会话的证据。",
                    "Figshare V1 的另一个归档映射是 DANDI 000166，未重复计入图谱；原 13 会话的独立动物数无法从已保存映射核实。",
                    "同场同步核团与 PRE/POST 保留共享动物。采用编码及来源级去重，无法排除未披露的原始动物重用。没有跨物种或跨数据集合并动物数。",
                    "Phase 10 和 Phase 11.1 的 bootstrap 单位是动物。历史 Phase 8/9 的置信区间按报告标明的动物重抽样，技术重复预测区间不是动物置信区间。"
                 ])
    write(ROOT/(OUT+"independence.json"),receipt)
    print(len(summaries),"cohorts,",len(checks),"source pairs,",sum(len(x) for x in assets.values()),"unique dataset-qualified recording IDs")


if __name__=="__main__":
    run()

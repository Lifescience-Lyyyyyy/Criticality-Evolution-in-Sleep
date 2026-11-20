#!/usr/bin/env python3
"""Render binned Allen CCFv3 label sections to outlines for the atlas."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from contourpy import contour_generator
from matplotlib.path import Path as PlotPath

from build_atlas_data import DEFAULT, write
from catalog import MOUSE

VOLUME = "https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_25.nrrd"
ONTOLOGY = "https://api.brain-map.org/api/v2/structure_graph_download/1.json"
EXPECTED = ("c620cbcc562183e4dcd40250d440130501781f74b41de35b1c1bdabace290c42",
            "b5f0d024d1df09ee18aef15093f55240f4a70b3924ee9689add612dadff80386")
BROAD = {
    "Isocortex": "#e5ecef", "OLF": "#ece8e4", "HPF": "#cfdde3", "CTXsp": "#e9e3de",
    "STR": "#dce5e5", "PAL": "#e5e8df", "TH": "#cedee2", "HY": "#dfdce3",
    "MB": "#ced8e1", "HB": "#d8dbe5", "CB": "#e5e4eb", "fiber tracts": "#f5f4f1", "VS": "#fafbfc",
}


def tree(entry):
    yield entry
    for child in entry["children"]:
        yield from tree(child)


def contours(mask):
    z = np.pad(mask, 1).astype("uint8")
    if not mask.any():
        return None
    segs = contour_generator(z=z, line_type="Separate").lines(.5)
    paths = []
    for point in segs:
        # At 50 μm sampling, drop redundant vertices and keep the volume-derived outline.
        p = PlotPath(point - .5, closed=True).cleaned(simplify=True)
        v = p.vertices[(p.codes == PlotPath.MOVETO) | (p.codes == PlotPath.LINETO)]
        if len(v) >= 3:
            coords = " ".join(f"{x:g},{y:g}" for x, y in v)
            paths.append(f"M{coords}Z")
    return "".join(paths) or None


def read_volume(path):
    raw = path.read_bytes()
    header, zipped = raw.split(b"\n\n", 1)
    for x in (b"sizes: 528 320 456", b"type: unsigned int", b"endian: little", b"encoding: gzip"):
        assert x in header
    # CCFv3 ReferenceSpace 9: AP, SI, LR, x-fastest in NRRD.
    volume = np.frombuffer(gzip.decompress(zipped), dtype="<u4").reshape((528, 320, 456), order="F")
    return np.ascontiguousarray(volume[::2, ::2, ::2])


def export(volume_file, ontology_file, out):
    for file, digest in zip((volume_file, ontology_file), EXPECTED, strict=True):
        assert hashlib.sha256(file.read_bytes()).hexdigest() == digest, file
    vol = read_volume(volume_file)
    data = json.loads(ontology_file.read_text())
    nodes = {it["acronym"]: it for it in tree(data["msg"][0])}
    ids = {r[0]: {c["id"] for c in tree(nodes[r[4]])} for r in MOUSE}
    groups = {name: {c["id"] for c in tree(nodes[name])} for name in BROAD}
    values, lookup = np.unique(vol, return_inverse=True)
    codes = lookup.reshape(vol.shape)
    del vol, lookup

    def mask(cut, selected):
        return np.isin(values, list(selected))[cut]

    results = {}
    for mode, axis in (("sagittal", 2), ("coronal", 0)):
        regions = {}
        for region, ancestors in ids.items():
            counts = mask(codes, ancestors).sum(axis=tuple(i for i in range(3) if i != axis))
            if mode == "sagittal":
                counts[len(counts) // 2 + 1:] = 0
            assert counts.max() > 0, (mode, region)
            regions[region] = int(counts.argmax())
        slices = {}
        for pos in sorted(set(regions.values())):
            cut = codes[:, :, pos].T if mode == "sagittal" else codes[pos, :, :]
            base = contours(cut != 0)
            labels = {}
            for name, sel in ids.items():
                m = mask(cut, sel)
                d = contours(m)
                if d:
                    yx = np.argwhere(m)
                    mid = yx.mean(axis=0)
                    y, x = yx[np.argmin(((yx - mid) ** 2).sum(axis=1))]
                    labels[name] = {"d": d, "center": [int(x), int(y)], "area": int(m.sum())}
            broad = [{"d": d, "fill": BROAD[name]} for name, sel in groups.items()
                     if (d := contours(mask(cut, sel)))]
            h, w = cut.shape
            slices[str(pos)] = {"w": w, "h": h, "outline": base, "broad": broad, "labels": labels}
        for region, pos in regions.items():
            assert region in slices[str(pos)]["labels"]
        results[mode] = {"positions": regions, "slices": slices}

    write(out, {
        "title": "Allen Mouse Brain CCFv3 · 2017",
        "source": VOLUME, "sourceInfo": "https://brain-map.org/support/documentation/api-allen-brain-connectivity-atlas",
        "ontology": ONTOLOGY, "paper": "https://doi.org/10.1016/j.cell.2020.04.007",
        "license": "https://alleninstitute.org/legal/terms-of-use",
        "sha256": {"volume": EXPECTED[0], "ontology": EXPECTED[1]},
        "voxelMm": .05, "species": "mouse", "ordering": ["AP", "SI", "LR"],
        "allenIds": {r[0]: {"acronym": r[4], "id": nodes[r[4]]["id"], "name": nodes[r[4]]["name"]} for r in MOUSE},
        "views": results,
    })


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("volume", type=Path)
    p.add_argument("ontology", type=Path)
    p.add_argument("--out", type=Path, default=DEFAULT / "anatomy.json")
    a = p.parse_args()
    export(a.volume, a.ontology, a.out)

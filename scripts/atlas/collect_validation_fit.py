"""Verify and unpack a Phase 11 MATLAB Actions result without persisting the signed URL."""
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from sys import argv
from time import sleep
from urllib.request import urlopen
from zipfile import ZipFile

from build_atlas_data import ROOT, write
from catalog import BASE, PIN
import json

run, stem, artifact, expected, url, *destination = argv[1:]
assert Path(stem).name == stem
for attempt in range(3):
    try:
        with urlopen(url, timeout=30) as stream:
            archive = stream.read()
        break
    except OSError:
        if attempt == 2:
            raise SystemExit("Unable to download Actions artifact after 3 attempts.")
        sleep(2)
assert "sha256:" + sha256(archive).hexdigest() == expected
with ZipFile(BytesIO(archive)) as z:
    assert z.namelist() == ["matlab_fits.json"] and not z.testzip()
    raw = z.read("matlab_fits.json")
doc = json.loads(raw)
assert str(doc["run_id"]) == run and doc["toolkit_commit"] == PIN and "R2024b" in doc["matlab"]
records = doc["records"]
assert len(records) == doc["input"]["records"]
for row in records:
    for f in [row["fit"], *([row["long_reach"]] if row["long_reach"] else [])]:
        assert f["exit_status"] == 0 and f["solver_flag"] > 0 and f["min_ar_root"] > 1 and f["d2"] >= 0
out = (ROOT / destination[0]) if destination else ROOT / "results/validation"
out.mkdir(parents=True,exist_ok=True)
(out / (stem + ".json")).write_bytes(raw)
write(out / (stem + ".receipt.json"), dict(
    action=f"{BASE}/actions/runs/{run}", artifact=artifact, artifact_sha256=expected,
    fit_sha256=sha256(raw).hexdigest(), commit=doc["project_commit"], matlab=doc["matlab"],
    count=doc["input"]["records"], input_sha256=doc["input"]["sha256"],
))
print("Verified run", run, "with", len(records), "rows.", flush=True)

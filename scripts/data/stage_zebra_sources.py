"""Transfer selected public GIN files from WSL to the mounted Colab Drive."""

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from hashlib import file_digest
from pathlib import Path
import json
import shlex
import shutil
import subprocess
from tempfile import TemporaryDirectory
from urllib.request import urlopen


paths = json.loads(Path("config/paths.json").read_text())
source = json.loads(Path(paths["source_manifest"]).read_text())["zebra_finch"]
settings = json.loads(Path("config/zebra_finch.json").read_text())
root = Path(paths["drive_root"])
gpu = ["/root/.local/bin/codex-gpu", "-s", "criticality-evolution"]

required = ["code/UtilsMatlab/birdParams.mat"]
for row in settings["recordings"]:
    base = f"dataset/sub-{row['bird']}/ses-{row['session']}"
    required += [f"{base}/ephys/{row['ap_meta']}", f"{base}/derived/{row['score']}"]
    required += [f"{base}/derived/{row['prefix']}{name}"
                 for name in ("cluster_info.tsv", "spike_times.npy", "spike_clusters.npy")]


def verify_release():
    sha = subprocess.check_output(["git", "ls-remote", source["repository"] + ".git", "HEAD"],
                                  text=True, timeout=30).split()[0]
    if sha != source["commit"]:
        raise ValueError(f"Expected GIN {source['commit']}, found {sha}")


def download(scratch, relpath):
    local = scratch / relpath
    local.parent.mkdir(parents=True, exist_ok=True)
    url = f"{source['repository']}/raw/master/{relpath}"
    with urlopen(url, timeout=75) as response, local.open("wb") as output:
        expected = response.headers.get("Content-Length")
        shutil.copyfileobj(response, output)
    size = local.stat().st_size
    if expected and size != int(expected):
        raise ValueError(f"GIN size mismatch for {relpath}")
    with local.open("rb") as reader:
        if local.suffix == ".npy" and reader.read(6) != b"\x93NUMPY":
            raise ValueError(f"Invalid NumPy data for {relpath}")
        reader.seek(0)
        digest = file_digest(reader, "sha256").hexdigest()
    print(f"Verified public GIN file, {size} bytes: {relpath}", flush=True)
    return relpath, {"source": url, "version": source["commit"], "bytes": size,
                     "sha256": digest, "scratch_path": str(local),
                     "drive_path": "datasets/zebra_finch/" + relpath}


verify_release()
with TemporaryDirectory(prefix="gin-stage-") as directory:
    scratch = Path(directory)
    with ThreadPoolExecutor(max_workers=4) as pool:
        assets = dict(pool.map(lambda name: download(scratch, name), required))
    verify_release()

    log = Path("logs")
    log.mkdir(exist_ok=True)
    folders = sorted({str((root / item["drive_path"]).parent) for item in assets.values()})
    mkdir = "mkdir -p " + " ".join(map(shlex.quote, folders))
    subprocess.run([*gpu, "run", mkdir, "--timeout", "120", "--output",
                    str(log / "gin_directories.ipynb")], check=True)
    for name, item in assets.items():
        destination = root / item["drive_path"]
        subprocess.run([*gpu, "upload", str(scratch / name), str(destination) + ".part"], check=True)

    record = {"version": source["commit"], "verified_in_wsl_at_utc": datetime.now(timezone.utc).isoformat(),
              "files": assets}
    transfer = scratch / "zebra-finch-transfer.json"
    transfer.write_text(json.dumps(record, indent=2) + "\n")
    receipt = root / "manifests" / transfer.name
    subprocess.run([*gpu, "upload", str(transfer), str(receipt)], check=True)

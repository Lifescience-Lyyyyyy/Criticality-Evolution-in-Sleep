"""Inventory official CRCNS sleep sources without fetching broadband recordings."""

from collections import defaultdict
from hashlib import file_digest
from io import RawIOBase
from pathlib import Path
from urllib.request import Request, urlopen
from zipfile import ZipFile, ZIP_DEFLATED
import argparse
import json
import re

HC23_SOURCE = "https://ndownloader.figshare.com/files/14834216"
HC23_SOURCE_BYTES = 493578348
HC23_MEMBER = re.compile(r"RAW_DATA/(Rat-\d+)/(Session\d+)/Day1/Overnight/(Sleep\d+)/TT(\d+)_(\d+)(_NonGC)?\.t")


class RemoteZip(RawIOBase):
    def __init__(self, url, size):
        self.url = url
        self.size = size
        self.position = 0
        self.start = 0
        self.cache = b""
        self.fetched_bytes = 0

    def seek(self, offset, whence=0):
        self.position = (0, self.position, self.size)[whence] + offset
        if self.position < 0:
            raise ValueError("Negative byte offset")
        return self.position

    def tell(self):
        return self.position

    def fetch(self, start, stop):
        end = stop - 1
        request = Request(self.url, headers={"Range": f"bytes={start}-{end}"})
        with urlopen(request, timeout=45) as response:
            expected = f"bytes {start}-{end}/{self.size}"
            if response.status != 206 or response.headers.get("Content-Range") != expected:
                raise ValueError(f"Invalid byte-range response for {expected}")
            data = response.read(stop - start + 1)
        if len(data) != stop - start:
            raise ValueError("Truncated byte range")
        self.fetched_bytes += len(data)
        return data

    def read(self, count=-1):
        end = self.size if count < 0 else min(self.size, self.position + count)
        length = end - self.position
        if length <= 0:
            return b""
        index = self.position - self.start
        if 0 <= index and index + length <= len(self.cache):
            data = self.cache[index:index + length]
        else:
            data = self.fetch(self.position, end)
        self.position = end
        return data

    def seekable(self):
        return True

    def readable(self):
        return True


def read_overnight(out):
    """Retrieve only sorted unit members; verify each source ZIP CRC."""
    receipt = out.with_suffix(".json")
    if out.exists() and receipt.exists():
        old = json.loads(receipt.read_text())
        with out.open("rb") as f:
            if file_digest(f, "sha256").hexdigest() == old["cache_sha256"]:
                return old

    fetched = RemoteZip(HC23_SOURCE, HC23_SOURCE_BYTES)
    progress = out.with_name(out.name + ".part")
    ranges = defaultdict(list)
    with ZipFile(fetched) as package, ZipFile(progress, "w", ZIP_DEFLATED) as stage:
        ordered = sorted(package.infolist(), key=lambda f: f.header_offset)
        for i, source in enumerate(ordered):
            match = HC23_MEMBER.fullmatch(source.filename)
            if match:
                stop = ordered[i + 1].header_offset if i + 1 < len(ordered) else package.start_dir
                ranges[match.group(1, 2)].append((source, stop))
        if len(ranges) != 29:
            raise ValueError(f"Unexpected source session inventory: {len(ranges)}")
        for subject_session, members in sorted(ranges.items()):
            start, end = members[0][0].header_offset, members[-1][1]
            if end - start > 8 * 2**20:
                raise ValueError(f"Unexpected overnight range in {subject_session}")
            fetched.start = start
            fetched.cache = fetched.fetch(start, end)
            for source, _ in members:
                stage.writestr(source.filename, package.read(source))
            print(*subject_session, len(members), end - start, flush=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    with progress.open("rb") as f:
        digest = file_digest(f, "sha256").hexdigest()
    result = {"source_file_id": 14834216, "source_zip_bytes": HC23_SOURCE_BYTES,
              "downloaded_bytes": fetched.fetched_bytes, "cache_sha256": digest,
              "sessions": len(ranges), "t_members": sum(len(value) for value in ranges.values())}
    progress.replace(out)
    receipt.write_text(json.dumps(result, indent=2) + "\n")
    return result


def audit_hc23(mat, source, config):
    from crcns_adapter import hc23_prepared
    from prepare_sleep_dandi import clock_bins

    width = config["bin_seconds"]
    records = []
    for r in hc23_prepared(mat, source):
        durations = {state: round(width * sum(z - a for a, z in clock_bins(ranges, width, config["min_segment_seconds"])), 8)
                     for state, ranges in r.sleep_intervals.items()}
        units = len(r.unit_ids)
        paired = units >= config["min_units"] and min(durations.values()) >= config["min_state_seconds"]
        item = {"dataset": r.source_dataset, "subject": r.subject_id, "session": r.session_id,
                "asset_id": r.asset_id, "region": r.region, "units": units, "state_seconds": durations,
                "recording_support_seconds": round(sum(b-a for a,b in r.recording_intervals),6),
                "total_support_segments": len(r.recording_intervals),
                "status": "PARTIAL" if paired else "BLOCKED", "eligible_comparisons": ["REM_minus_NREM"] if paired else [],
                "reason": "Author-validated REM/NREM; overnight Wake unavailable" if paired
                          else "Insufficient verified stable units or REM/NREM exposure",
                "quality": r.quality_metrics, "clock": r.timebase, "provenance": r.provenance}
        print(r.subject_id, r.session_id, units, durations, item["status"], flush=True)
        records.append(item)
    return records


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hc23-overnight", required=True, type=Path)
    p.add_argument("--hc23-mat", type=Path)
    p.add_argument("--config", type=Path, default=Path("config/sleep_dandi.json"))
    p.add_argument("--out", type=Path)
    args = p.parse_args()
    args.hc23_overnight.parent.mkdir(parents=True, exist_ok=True)
    source = read_overnight(args.hc23_overnight)
    print(source)
    if args.hc23_mat:
        if not args.out:
            p.error("--out is required when auditing the MAT")
        records = audit_hc23(args.hc23_mat, args.hc23_overnight, json.loads(args.config.read_text()))
        archive = {"source": source, "records": records}
        temporary = args.out.with_suffix(".part")
        temporary.write_text(json.dumps(archive, allow_nan=False, indent=2) + "\n")
        temporary.replace(args.out)

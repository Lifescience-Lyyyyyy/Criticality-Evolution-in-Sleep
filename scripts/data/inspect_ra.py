"""Inventory Elmaleh et al.'s public OSF release without downloading recordings."""

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
import struct
from urllib.parse import urlparse
from urllib.request import Request, urlopen
import zlib


NODE = "https://api.osf.io/v2/nodes/8wfz6/"
MAX_SLICE = 8192


def api(url):
    if urlparse(url).hostname != "api.osf.io":
        raise ValueError(f"Unexpected OSF metadata URL: {url}")
    with urlopen(Request(url, headers={"Accept": "application/json"}), timeout=30) as response:
        return json.load(response)


def list_rows(url):
    while url:
        result = api(url)
        yield from result["data"]
        url = result["links"]["next"]


def items(row):
    attr = row["attributes"]
    if attr["kind"] == "folder":
        child_url = row["relationships"]["files"]["links"]["related"]["href"]
        for child in list_rows(child_url):
            yield from items(child)
    else:
        yield {"id": row["id"], "metadata": row["links"]["self"], "download": row["links"]["download"],
               "path": attr["materialized_path"], "bytes": attr["size"], "version": attr["current_version"],
               "modified": attr["date_modified"], "osf_sha256": attr["extra"]["hashes"]["sha256"]}


def read_range(row, offset, limit):
    end = min(offset + limit, row["bytes"]) - 1
    if not 0 <= offset <= end or not 0 < limit <= MAX_SLICE:
        raise ValueError(f"Invalid bounded request for {row['path']}: {offset}, {limit}")
    header = f"bytes={offset}-{end}"
    expected = f"bytes {offset}-{end}/{row['bytes']}"
    request = Request(row["download"], headers={"Range": header, "Accept-Encoding": "identity"})
    with urlopen(request, timeout=30) as response:
        if response.status != 206 or response.headers.get("Content-Range") != expected:
            raise ValueError(f"{row['path']}: expected HTTP 206 {expected}, received "
                             f"{response.status} {response.headers.get('Content-Range')}")
        raw = response.read(end - offset + 2)
    if len(raw) != end - offset + 1:
        raise ValueError(f"{row['path']}: truncated range {header}")
    return raw


def mat_elements(raw, offset=0):
    while offset + 8 <= len(raw):
        tag, length = struct.unpack_from("<II", raw, offset)
        if tag >> 16:
            yield tag & 65535, raw[offset + 4:offset + 4 + (tag >> 16)]
            offset += 8
        else:
            yield tag, raw[offset + 8:offset + 8 + length]
            offset += 8 + (length + 7) // 8 * 8


def array_header(raw):
    elements = mat_elements(raw)
    _, flags = next(elements)
    _, shape = next(elements)
    _, name = next(elements)
    return {"name": name.decode("ascii"), "shape": [n for n, in struct.iter_unpack("<i", shape)],
            "matlab_class": flags[0]}, elements


def mat_schema(row):
    initial = read_range(row, 0, MAX_SLICE)
    if not initial.startswith(b"MATLAB 5.0 MAT-file") or initial[126:128] != b"IM":
        raise ValueError(f"{row['path']}: expected little-endian MATLAB Level 5")
    variables = []
    hemisphere_indices = {}
    offset = 128
    while offset < row["bytes"]:
        if len(variables) >= 128:
            raise ValueError(f"{row['path']}: MAT metadata scan exceeded 128 variables")
        block = initial[offset:] if offset + 4096 <= len(initial) else read_range(row, offset, 4096)
        tag, length = struct.unpack_from("<II", block)
        if tag == 15:
            # Decode only the header and a small numeric prefix. Skip the compressed payload.
            inner = zlib.decompressobj().decompress(block[8:8 + length], 32768)
            _, matrix = next(mat_elements(inner))
            offset += 8 + length
        elif tag == 14:
            matrix = block[8:8 + length]
            offset += 8 + (length + 7) // 8 * 8
        else:
            raise ValueError(f"{row['path']}: unexpected MATLAB element type {tag} at byte {offset}")
        info, elements = array_header(matrix)
        if info["matlab_class"] == 2:
            _, field_length = next(elements)
            width, = struct.unpack("<i", field_length)
            _, names = next(elements)
            info["fields"] = [names[i:i + width].rstrip(b"\0").decode("ascii")
                              for i in range(0, len(names), width)]
            if info["fields"][0] == "spikeTimes":
                _, first_field = next(elements)
                first, values = array_header(first_field)
                storage_type, samples = next(values)
                dtype = {9: "d", 7: "f"}[storage_type]
                info["first_cluster_spike_vector_shape"] = first["shape"]
                info["first_cluster_first_four_times"] = list(struct.unpack("<4" + dtype, samples[:4 * struct.calcsize(dtype)]))
        if info["name"] in ("indL", "indR"):
            storage_type, data = next(elements)
            dtype = {1: "b", 2: "B", 3: "h", 4: "H", 5: "i", 6: "I", 9: "d"}[storage_type]
            indices = [n for n, in struct.iter_unpack("<" + dtype, data)]
            if any(n != int(n) for n in indices):
                raise ValueError(f"{row['path']}: hemisphere array has noninteger indices")
            hemisphere_indices[info["name"]] = indices
        variables.append(info)
    if offset != row["bytes"]:
        raise ValueError(f"{row['path']}: MATLAB top-level elements did not reach recorded size")
    rows = {v["name"]: v for v in variables}
    units = rows["clusters"]["shape"][1]
    left, right = hemisphere_indices["indL"], hemisphere_indices["indR"]
    if sorted(left + right) != list(range(1, units + 1)):
        raise ValueError(f"{row['path']}: indL and indR do not partition clusters")
    return variables


def audio_schema(row):
    raw = read_range(row, 0, 4096)
    if raw[:4] != b"RIFF" or raw[8:12] != b"WAVE":
        raise ValueError(f"{row['path']}: expected RIFF WAVE")
    chunks = []
    offset = 12
    fmt = None
    data_bytes = None
    while offset + 8 <= len(raw):
        label, length = struct.unpack_from("<4sI", raw, offset)
        chunks.append(label.decode("ascii"))
        if label == b"fmt ":
            fmt = struct.unpack_from("<HHIIHH", raw, offset + 8)
        elif label == b"data":
            data_bytes = length
        offset += 8 + (length + 1) // 2 * 2
    if fmt is None or data_bytes is None:
        raise ValueError(f"{row['path']}: format and data chunks were not in the inspected prefix")
    code, channels, hz, byte_rate, _, bits = fmt
    return {"chunks": chunks, "format": code, "channels": channels, "sample_rate_hz": hz, "bits": bits,
            "audio_duration_s": data_bytes / byte_rate, "bytes_after_inspected_chunks": row["bytes"] - offset}


def describe(row):
    if row["path"].endswith(".mat"):
        return {**row, "mat_variables": mat_schema(row), "hemisphere_indices_partition_clusters": True}
    if row["path"].endswith(".wav"):
        return {**row, "audio_header": audio_schema(row)}
    return row


def inspect():
    attr = api(NODE)["data"]["attributes"]
    if not attr["public"]:
        raise ValueError("OSF 8wfz6 is not publicly readable")
    providers = list(list_rows(NODE + "files/"))
    if [provider["attributes"]["name"] for provider in providers] != ["osfstorage"]:
        raise ValueError("OSF storage providers changed. Recheck the inventory.")
    roots = list(list_rows(NODE + "files/osfstorage/"))
    with ThreadPoolExecutor(max_workers=5) as pool:
        files = sorted((item for group in pool.map(items, roots) for item in group), key=lambda r: r["path"])
        evidence = list(pool.map(describe, files))
    return {"schema_version": 1, "audited_utc": datetime.now(timezone.utc).isoformat(),
            "doi": "10.17605/OSF.IO/8WFZ6", "project_api": NODE, "public": True,
            "license": attr["node_license"], "metadata_counts": {
                kind: len(list(list_rows(NODE + kind + "/"))) for kind in ("children", "linked_nodes", "wikis")},
            "sampling": "OSF file metadata plus bounded HTTP 206 reads, no complete recordings",
            "files": evidence}


if __name__ == "__main__":
    print(json.dumps(inspect(), indent=2))

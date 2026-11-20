"""Read selected Zenodo members with central-directory size and CRC verification."""
import argparse
from pathlib import Path, PurePosixPath
import json
from struct import unpack_from
from time import sleep
from urllib.error import HTTPError, URLError
import zlib

from audit_crcns_sleep import RemoteZip


def member(remote,item):
    name=item["name"]
    start=item["header_offset"]
    stop=min(remote.size,start+30+len(name.encode())+4096+item["compress_size"])
    for n in range(4):
        try:
            raw=remote.fetch(start,stop)
            break
        except (TimeoutError,URLError) as e:
            if n==3 or isinstance(e,HTTPError) and e.code not in (429,500,502,503,504):
                raise
            retry=e.headers.get("Retry-After") if isinstance(e,HTTPError) else None
            sleep(min(int(retry) if retry and retry.isdigit() else 15*(n+1),120))
    signature,_,_,method=unpack_from("<4sHHH",raw)
    label,extra=unpack_from("<HH",raw,26)
    offset=30+label+extra
    data=raw[offset:offset+item["compress_size"]]
    if signature!=b"PK\x03\x04" or raw[30:30+label].decode()!=name or len(data)!=item["compress_size"]:
        raise ValueError(f"Local ZIP header disagrees with central directory: {name}")
    if method not in (0,8):
        raise ValueError(f"Unsupported ZIP compression method {method}")
    content=zlib.decompress(data,-15) if method==8 else data
    if len(content)!=item["file_size"] or zlib.crc32(content)!=item["crc"]:
        raise ValueError(f"ZIP CRC or size mismatch for {name}")
    return content


def download(metadata,inventory,dest,names):
    source,=metadata["files"]
    entries={row["name"]:row for row in inventory}
    remote=RemoteZip(source["links"]["self"],source["size"])
    for name in names:
        rel=PurePosixPath(name).relative_to("lozano_hong_jin_2026")
        if ".." in rel.parts or not rel.parts:
            raise ValueError(f"Unsafe source member {name}")
        path=dest/rel
        item=entries[name]
        if path.exists() and path.stat().st_size==item["file_size"] and zlib.crc32(path.read_bytes())==item["crc"]:
            continue
        content=member(remote,item)
        path.parent.mkdir(parents=True,exist_ok=True)
        partial=path.with_suffix(path.suffix+".part")
        partial.write_bytes(content)
        partial.replace(path)
        print(name,len(content),flush=True)
    print("CRC-verified range bytes:",remote.fetched_bytes,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-meta",type=Path,required=True)
    p.add_argument("--inventory",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("members",nargs="+")
    a=p.parse_args()
    download(json.loads(a.source_meta.read_text()),json.loads(a.inventory.read_text()),a.out,a.members)

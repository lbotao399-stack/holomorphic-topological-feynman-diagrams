"""Store every exact relation in deterministic gzip, without dropping rows."""
from hashlib import sha256
import gzip
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def compact(path):
    data = json.loads(path.read_text())
    if "rows" not in data:
        return
    rows = data.pop("rows")
    raw = (json.dumps(rows, separators=(",", ":")) + "\n").encode()
    packed = gzip.compress(raw, compresslevel=9, mtime=0)
    target = path.with_suffix(".rows.json.gz")
    target.write_bytes(packed)
    if gzip.decompress(target.read_bytes()) != raw:
        raise ArithmeticError("relation compression changed the exact row data")
    data["relation_file"] = target.name
    data["relation_count"] = len(rows)
    data["relation_uncompressed_sha256"] = sha256(raw).hexdigest()
    data["relation_compressed_sha256"] = sha256(packed).hexdigest()
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    temporary.replace(path)
    print(path.name, len(rows), "complete rows", len(raw), "->", len(packed), "bytes")


if __name__ == "__main__":
    for path in sorted((HERE / "tables").glob("*.json")):
        compact(path)

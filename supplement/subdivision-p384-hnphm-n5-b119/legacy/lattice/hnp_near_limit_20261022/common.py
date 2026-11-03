"""Small immutable-record helpers shared by generation and execution."""
import hashlib
import json
import os
from pathlib import Path
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    with Path(path).open("x") as handle:
        json.dump(value, handle, sort_keys=True, indent=2)
        handle.write("\n")


def atomic(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".writing")
    with temporary.open("w") as handle:
        json.dump(value, handle, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def utc():
    return datetime.now(timezone.utc).isoformat()


def uniform_below(bound, material):
    bits = bound.bit_length()
    counter = 0
    while True:
        data = b"".join(bytes.fromhex(digest([material, counter, block]))
                        for block in range((bits + 255) // 256))
        value = int.from_bytes(data, "big") & ((1 << bits) - 1)
        if value < bound:
            return value
        counter += 1

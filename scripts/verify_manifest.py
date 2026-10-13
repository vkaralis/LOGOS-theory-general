"""Verify release SHA-256 file hashes from any working directory."""
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'MANIFEST_SHA256.json').read_text())
for name,expected in manifest.items():
    file=ROOT/name
    if not file.is_file() or hashlib.sha256(file.read_bytes()).hexdigest()!=expected:
        raise SystemExit(f'Integrity failure: {name}')
print(f'PASS: {len(manifest)} release file hashes')

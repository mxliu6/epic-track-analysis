#!/usr/bin/env python3
"""Validate ROOT integrity and required ePIC tracking collections before a campaign run."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import uproot


REQUIRED = (
    "MCParticles/MCParticles.PDG",
    "MCParticles/MCParticles.generatorStatus",
    "MCParticles/MCParticles.charge",
    "ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.x",
    "_ReconstructedChargedParticleAssociations_rec/_ReconstructedChargedParticleAssociations_rec.index",
    "_ReconstructedChargedParticleAssociations_sim/_ReconstructedChargedParticleAssociations_sim.index",
)


def header_end(path: Path) -> int | None:
    with path.open("rb") as stream:
        header = stream.read(16)
    if len(header) < 16 or header[:4] != b"root":
        return None
    version = struct.unpack(">i", header[4:8])[0]
    if version >= 1_000_000:
        return None  # Large-file ROOT headers use a different field layout.
    return struct.unpack(">i", header[12:16])[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path)
    args = parser.parse_args()
    failures = 0
    for path in args.files:
        expected = header_end(path)
        actual = path.stat().st_size
        if expected is not None and expected > actual:
            failures += 1
            print(f"BAD  {path.name}: truncated by {expected - actual:,} bytes (header expects {expected:,})")
            continue
        try:
            with uproot.open(path) as root_file:
                events = root_file["events"]
                names = set(events.keys(recursive=True))
                missing = [name for name in REQUIRED if name not in names]
                if missing:
                    failures += 1
                    print(f"BAD  {path.name}: missing {len(missing)} required branches")
                else:
                    print(f"OK   {path.name}: {events.num_entries:,} events")
        except Exception as exc:
            failures += 1
            print(f"BAD  {path.name}: {type(exc).__name__}: {str(exc).splitlines()[0]}")
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()


#!/usr/bin/env python3
"""Inspect an ePIC EDM4eic ROOT file without requiring CERN ROOT."""

from __future__ import annotations

import argparse
from pathlib import Path

import uproot


INTERESTING = (
    "MCParticles",
    "ReconstructedChargedParticles",
    "ReconstructedChargedParticleAssociations",
    "_ReconstructedChargedParticleAssociations",
    "ReconstructedChargedParticleLinks",
    "_ReconstructedChargedParticleLinks",
    "podio_metadata",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path)
    parser.add_argument("--all", action="store_true", help="show every event branch")
    args = parser.parse_args()

    with uproot.open(args.file) as root_file:
        print(f"File: {args.file}")
        print("Top-level objects:")
        for key, classname in root_file.classnames().items():
            print(f"  {key:<40} {classname}")

        if "events" not in root_file:
            print("\nNo 'events' object was found. This may be generator-level or DD4hep output,")
            print("rather than EICrecon EDM4eic reconstruction output.")
            return

        events = root_file["events"]
        print(f"\nEvents: {events.num_entries:,}")
        print("Relevant event branches:")
        keys = events.keys(recursive=True)
        selected = keys if args.all else [k for k in keys if k.startswith(INTERESTING)]
        for key in selected:
            branch = events[key]
            typename = getattr(branch, "typename", type(branch).__name__)
            print(f"  {key:<72} {typename}")

        names = set(keys)
        required_groups = {
            "truth particles": ["MCParticles/MCParticles.PDG", "MCParticles/MCParticles.momentum.x"],
            "charged reconstructed particles": [
                "ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.x",
                "ReconstructedChargedWithoutPIDParticles/ReconstructedChargedWithoutPIDParticles.momentum.x",
            ],
            "truth associations": [
                "_ReconstructedChargedParticleAssociations_rec/_ReconstructedChargedParticleAssociations_rec.index",
                "_ReconstructedChargedWithoutPIDParticleAssociations_rec/_ReconstructedChargedWithoutPIDParticleAssociations_rec.index",
            ],
        }
        print("\nAnalysis readiness:")
        for label, alternatives in required_groups.items():
            found = next((name for name in alternatives if name in names), None)
            print(f"  {'yes' if found else 'NO ':<4} {label}" + (f" ({found})" if found else ""))


if __name__ == "__main__":
    main()

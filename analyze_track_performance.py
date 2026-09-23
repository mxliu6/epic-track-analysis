#!/usr/bin/env python3
"""Measure charged-track efficiency, purity, duplicates, and PID response in EDM4eic files."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from collections import defaultdict
from pathlib import Path

# Keep plotting caches local and writable on managed macOS environments.
_CACHE = Path(__file__).resolve().parent / ".cache"
(_CACHE / "matplotlib").mkdir(parents=True, exist_ok=True)
(_CACHE / "fontconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE))
os.environ.setdefault("MPLBACKEND", "Agg")

import awkward as ak
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import uproot


SPECIES = {11: "electron_or_positron", 13: "muon", 211: "pion", 321: "kaon", 2212: "proton"}
ELECTRON_CHARGES = {11: "electron", -11: "positron"}
RECO_COLLECTIONS = (
    ("ReconstructedChargedParticles", "ReconstructedChargedParticleAssociations"),
    ("ReconstructedChargedWithoutPIDParticles", "ReconstructedChargedWithoutPIDParticleAssociations"),
)


def branch_name(names: set[str], alternatives: list[str], label: str, required: bool = True) -> str | None:
    for name in alternatives:
        if name in names:
            return name
    if required:
        raise RuntimeError(f"Cannot find {label}; tried: {', '.join(alternatives)}")
    return None


def discover(path: str) -> dict[str, str | None]:
    with uproot.open(path) as root_file:
        if "events" not in root_file:
            raise RuntimeError(f"{path} has no 'events' tree/RNTuple")
        names = set(root_file["events"].keys(recursive=True))

    reco = assoc = None
    for reco_candidate, assoc_candidate in RECO_COLLECTIONS:
        if f"{reco_candidate}/{reco_candidate}.momentum.x" in names:
            reco, assoc = reco_candidate, assoc_candidate
            break
    if reco is None:
        raise RuntimeError("No supported reconstructed charged-particle collection was found")

    result: dict[str, str | None] = {
        "mc_pdg": branch_name(names, ["MCParticles/MCParticles.PDG"], "MC PDG"),
        "mc_status": branch_name(names, ["MCParticles/MCParticles.generatorStatus"], "MC generator status"),
        "mc_charge": branch_name(names, ["MCParticles/MCParticles.charge"], "MC charge"),
        "mc_px": branch_name(names, ["MCParticles/MCParticles.momentum.x"], "MC px"),
        "mc_py": branch_name(names, ["MCParticles/MCParticles.momentum.y"], "MC py"),
        "mc_pz": branch_name(names, ["MCParticles/MCParticles.momentum.z"], "MC pz"),
        "reco_px": branch_name(names, [f"{reco}/{reco}.momentum.x"], "reco px"),
        "reco_py": branch_name(names, [f"{reco}/{reco}.momentum.y"], "reco py"),
        "reco_pz": branch_name(names, [f"{reco}/{reco}.momentum.z"], "reco pz"),
        "reco_pdg": branch_name(names, [f"{reco}/{reco}.PDG"], "reco PDG", required=False),
        "assoc_rec": branch_name(
            names,
            [
                f"_{assoc}_rec/_{assoc}_rec.index",
                f"_{assoc}_recID/_{assoc}_recID.index",
            ],
            "association reco index",
        ),
        "assoc_sim": branch_name(
            names,
            [
                f"_{assoc}_sim/_{assoc}_sim.index",
                f"_{assoc}_simID/_{assoc}_simID.index",
            ],
            "association sim index",
        ),
        "assoc_weight": branch_name(names, [f"{assoc}/{assoc}.weight"], "association weight", required=False),
    }
    result["reco_collection"] = reco
    result["assoc_collection"] = assoc
    return result


def eta(px: np.ndarray, py: np.ndarray, pz: np.ndarray) -> np.ndarray:
    pt = np.hypot(px, py)
    return np.arcsinh(np.divide(pz, pt, out=np.full_like(pz, np.inf, dtype=float), where=pt > 0))


def binomial(numerator: np.ndarray, denominator: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    value = np.divide(numerator, denominator, out=np.full_like(numerator, np.nan, dtype=float), where=denominator > 0)
    error = np.sqrt(np.divide(value * (1.0 - value), denominator, out=np.zeros_like(value), where=denominator > 0))
    return value, error


def plot_curves(output: Path, bins: np.ndarray, curves: dict[str, tuple[np.ndarray, np.ndarray]], ylabel: str) -> None:
    centers = 0.5 * (bins[:-1] + bins[1:])
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    for label, (num, den) in curves.items():
        value, error = binomial(num, den)
        ax.errorbar(centers, value, yerr=error, marker="o", ms=3, capsize=2, label=label)
    ax.set(xlabel=r"$p_T$ [GeV]", ylabel=ylabel, ylim=(0, 1.05))
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=9)
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def resolution_statistics(values: list[list[float]], min_entries: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    count = np.asarray([len(v) for v in values], dtype=int)
    median = np.full(len(values), np.nan)
    sigma68 = np.full(len(values), np.nan)
    mean = np.full(len(values), np.nan)
    rms = np.full(len(values), np.nan)
    for i, entries in enumerate(values):
        if len(entries) < min_entries:
            continue
        data = np.asarray(entries)
        q16, q50, q84 = np.quantile(data, [0.16, 0.50, 0.84])
        median[i] = q50
        sigma68[i] = 0.5 * (q84 - q16)
        mean[i] = np.mean(data)
        rms[i] = np.std(data)
    return count, median, sigma68, mean, rms


def plot_resolution(
    output: Path,
    bins: np.ndarray,
    residuals: dict[str, list[list[float]]],
    min_entries: int,
    residual_label: str,
) -> None:
    centers = 0.5 * (bins[:-1] + bins[1:])
    fig, (bias_ax, width_ax) = plt.subplots(2, 1, figsize=(7.2, 7.2), sharex=True)
    for species, values in residuals.items():
        _, median, sigma68, _, _ = resolution_statistics(values, min_entries)
        bias_ax.plot(centers, median, marker="o", ms=3, label=species)
        width_ax.plot(centers, sigma68, marker="o", ms=3, label=species)
    bias_ax.axhline(0, color="black", lw=0.8)
    bias_ax.set_ylabel(f"Median {residual_label}")
    width_ax.set(xlabel=r"True $p_T$ [GeV]", ylabel=r"Resolution $\sigma_{68}$")
    for ax in (bias_ax, width_ax):
        ax.grid(alpha=0.25)
    bias_ax.legend(ncol=2, fontsize=9)
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", help="EICrecon EDM4eic ROOT file(s)")
    parser.add_argument("-o", "--output", type=Path, default=Path("track_performance"))
    parser.add_argument("--pt-bins", default="0,0.2,0.5,1,2,3,5,7.5,10,15,20,30,50")
    parser.add_argument("--eta-min", type=float, default=-3.5)
    parser.add_argument("--eta-max", type=float, default=3.5)
    parser.add_argument("--min-weight", type=float, default=0.0)
    parser.add_argument("--min-resolution-entries", type=int, default=5)
    parser.add_argument("--step-size", default="100 MB")
    args = parser.parse_args()

    bins = np.asarray([float(x) for x in args.pt_bins.split(",")])
    if len(bins) < 2 or np.any(np.diff(bins) <= 0):
        raise SystemExit("--pt-bins must be a strictly increasing comma-separated list")
    args.output.mkdir(parents=True, exist_ok=True)

    schema = discover(args.files[0])
    expressions = [v for k, v in schema.items() if k not in {"reco_collection", "assoc_collection"} and v]
    counts = {
        name: {"truth": np.zeros(len(bins) - 1, int), "found": np.zeros(len(bins) - 1, int)}
        for name in SPECIES.values()
    }
    electron_charge_counts = {
        name: {"truth": np.zeros(len(bins) - 1, int), "found": np.zeros(len(bins) - 1, int)}
        for name in ELECTRON_CHARGES.values()
    }
    purity_den = np.zeros(len(bins) - 1, int)
    purity_num = np.zeros(len(bins) - 1, int)
    pt_response = np.zeros((len(bins) - 1, len(bins) - 1), int)
    pid_den = {name: np.zeros(len(bins) - 1, int) for name in SPECIES.values()}
    pid_num = {name: np.zeros(len(bins) - 1, int) for name in SPECIES.values()}
    confusion: defaultdict[tuple[int, int], int] = defaultdict(int)
    pt_residuals = {
        name: [[] for _ in range(len(bins) - 1)] for name in SPECIES.values()
    }
    momentum_residuals = {
        name: [[] for _ in range(len(bins) - 1)] for name in SPECIES.values()
    }
    electron_charge_pt_residuals = {
        name: [[] for _ in range(len(bins) - 1)] for name in ELECTRON_CHARGES.values()
    }
    electron_charge_momentum_residuals = {
        name: [[] for _ in range(len(bins) - 1)] for name in ELECTRON_CHARGES.values()
    }
    totals = defaultdict(int)

    sources = [f"{path}:events" for path in args.files]
    for arrays in uproot.iterate(sources, expressions=expressions, step_size=args.step_size, library="ak"):
        nevents = len(arrays[schema["mc_pdg"]])
        totals["events"] += nevents
        for iev in range(nevents):
            def event(key: str, dtype: type = float) -> np.ndarray:
                return np.asarray(ak.to_list(arrays[schema[key]][iev]), dtype=dtype)

            mc_pdg = event("mc_pdg", int)
            mc_status = event("mc_status", int)
            mc_charge = event("mc_charge")
            mc_px, mc_py, mc_pz = event("mc_px"), event("mc_py"), event("mc_pz")
            mc_pt = np.hypot(mc_px, mc_py)
            mc_eta = eta(mc_px, mc_py, mc_pz)
            eligible = (mc_status == 1) & (np.abs(mc_charge) > 0) & (mc_eta >= args.eta_min) & (mc_eta <= args.eta_max)

            reco_px, reco_py, reco_pz = event("reco_px"), event("reco_py"), event("reco_pz")
            reco_pt = np.hypot(reco_px, reco_py)
            reco_eta = eta(reco_px, reco_py, reco_pz)
            reco_accept = (reco_eta >= args.eta_min) & (reco_eta <= args.eta_max)
            rec_idx, sim_idx = event("assoc_rec", int), event("assoc_sim", int)
            valid = (rec_idx >= 0) & (rec_idx < len(reco_pt)) & (sim_idx >= 0) & (sim_idx < len(mc_pdg))
            weights = event("assoc_weight") if schema["assoc_weight"] else np.ones(len(rec_idx))
            if schema["assoc_weight"]:
                valid &= weights >= args.min_weight
            pairs = {(int(r), int(s)) for r, s in zip(rec_idx[valid], sim_idx[valid])}
            best_truth: dict[int, int] = {}
            best_weight: dict[int, float] = {}
            for r, s, w, ok in zip(rec_idx, sim_idx, weights, valid):
                if ok and eligible[s] and (int(r) not in best_weight or w > best_weight[int(r)]):
                    best_truth[int(r)], best_weight[int(r)] = int(s), float(w)
            matched_rec = {r for r, s in pairs if eligible[s]}
            matched_sim = {s for r, s in pairs if reco_accept[r]}

            totals["truth_eligible"] += int(np.count_nonzero(eligible))
            totals["reco_accepted"] += int(np.count_nonzero(reco_accept))
            totals["reco_matched"] += sum(bool(reco_accept[r]) for r in matched_rec)
            sim_multiplicity = defaultdict(set)
            for r, s in pairs:
                if eligible[s] and reco_accept[r]:
                    sim_multiplicity[s].add(r)
            totals["duplicate_excess"] += sum(max(0, len(rs) - 1) for rs in sim_multiplicity.values())
            totals["matched_truth_eligible"] += len(sim_multiplicity)

            purity_den += np.histogram(reco_pt[reco_accept], bins=bins)[0]
            matched_mask = np.zeros(len(reco_pt), dtype=bool)
            if matched_rec:
                matched_mask[list(matched_rec)] = True
            purity_num += np.histogram(reco_pt[reco_accept & matched_mask], bins=bins)[0]
            response_rec = [r for r in best_truth if reco_accept[r]]
            if response_rec:
                response_sim = [best_truth[r] for r in response_rec]
                pt_response += np.histogram2d(mc_pt[response_sim], reco_pt[response_rec], bins=(bins, bins))[0].astype(int)
                mc_p = np.sqrt(mc_px**2 + mc_py**2 + mc_pz**2)
                reco_p = np.sqrt(reco_px**2 + reco_py**2 + reco_pz**2)
                for r, s in zip(response_rec, response_sim):
                    species = SPECIES.get(abs(int(mc_pdg[s])))
                    ibin = int(np.searchsorted(bins, mc_pt[s], side="right") - 1)
                    if species is None or not 0 <= ibin < len(bins) - 1 or mc_pt[s] <= 0 or mc_p[s] <= 0:
                        continue
                    pt_residuals[species][ibin].append(float((reco_pt[r] - mc_pt[s]) / mc_pt[s]))
                    momentum_residuals[species][ibin].append(float((reco_p[r] - mc_p[s]) / mc_p[s]))
                    charge_species = ELECTRON_CHARGES.get(int(mc_pdg[s]))
                    if charge_species:
                        electron_charge_pt_residuals[charge_species][ibin].append(
                            float((reco_pt[r] - mc_pt[s]) / mc_pt[s])
                        )
                        electron_charge_momentum_residuals[charge_species][ibin].append(
                            float((reco_p[r] - mc_p[s]) / mc_p[s])
                        )

            for pdg, name in SPECIES.items():
                truth_mask = eligible & (np.abs(mc_pdg) == pdg)
                found_mask = truth_mask & np.asarray([i in matched_sim for i in range(len(mc_pdg))])
                counts[name]["truth"] += np.histogram(mc_pt[truth_mask], bins=bins)[0]
                counts[name]["found"] += np.histogram(mc_pt[found_mask], bins=bins)[0]
            for pdg, name in ELECTRON_CHARGES.items():
                truth_mask = eligible & (mc_pdg == pdg)
                found_mask = truth_mask & np.asarray([i in matched_sim for i in range(len(mc_pdg))])
                electron_charge_counts[name]["truth"] += np.histogram(mc_pt[truth_mask], bins=bins)[0]
                electron_charge_counts[name]["found"] += np.histogram(mc_pt[found_mask], bins=bins)[0]

            if schema["reco_pdg"]:
                reco_pdg = event("reco_pdg", int)
                for r in np.flatnonzero(reco_accept):
                    reco_abs = abs(int(reco_pdg[r]))
                    if reco_abs in SPECIES:
                        name = SPECIES[reco_abs]
                        pid_den[name] += np.histogram([reco_pt[r]], bins=bins)[0]
                        if r in best_truth:
                            true_abs = abs(int(mc_pdg[best_truth[r]]))
                            confusion[(true_abs, reco_abs)] += 1
                            if true_abs == reco_abs:
                                pid_num[name] += np.histogram([reco_pt[r]], bins=bins)[0]

    plot_curves(
        args.output / "tracking_efficiency.png",
        bins,
        {name: (data["found"], data["truth"]) for name, data in counts.items()},
        "Tracking efficiency",
    )
    plot_curves(
        args.output / "electron_positron_tracking_efficiency.png",
        bins,
        {name: (data["found"], data["truth"]) for name, data in electron_charge_counts.items()},
        "Tracking efficiency",
    )
    plot_curves(args.output / "track_purity.png", bins, {"all charged": (purity_num, purity_den)}, "Matched-track purity")
    if schema["reco_pdg"]:
        plot_curves(args.output / "pid_purity.png", bins, {name: (pid_num[name], pid_den[name]) for name in SPECIES.values()}, "PID purity")
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    positive = pt_response[pt_response > 0]
    norm = LogNorm(vmin=1, vmax=max(1, int(positive.max()))) if len(positive) else None
    mesh = ax.pcolormesh(bins, bins, pt_response.T, shading="auto", norm=norm)
    ax.plot([bins[0], bins[-1]], [bins[0], bins[-1]], "w--", lw=1)
    ax.set(xlabel=r"True $p_T$ [GeV]", ylabel=r"Reconstructed $p_T$ [GeV]", xlim=(bins[0], bins[-1]), ylim=(bins[0], bins[-1]))
    fig.colorbar(mesh, ax=ax, label="Matched tracks")
    fig.tight_layout()
    fig.savefig(args.output / "pt_response.png", dpi=160)
    plt.close(fig)
    plot_resolution(
        args.output / "pt_resolution_vs_pt.png",
        bins,
        pt_residuals,
        args.min_resolution_entries,
        r"$\Delta p_T/p_T$",
    )
    plot_resolution(
        args.output / "momentum_resolution_vs_pt.png",
        bins,
        momentum_residuals,
        args.min_resolution_entries,
        r"$\Delta p/p$",
    )
    plot_resolution(
        args.output / "electron_positron_pt_resolution_vs_pt.png",
        bins,
        electron_charge_pt_residuals,
        args.min_resolution_entries,
        r"$\Delta p_T/p_T$",
    )
    plot_resolution(
        args.output / "electron_positron_momentum_resolution_vs_pt.png",
        bins,
        electron_charge_momentum_residuals,
        args.min_resolution_entries,
        r"$\Delta p/p$",
    )

    with (args.output / "binned_metrics.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["metric", "species", "pt_low", "pt_high", "numerator", "denominator", "value", "stat_error"])
        rows = [("tracking_efficiency", name, data["found"], data["truth"]) for name, data in counts.items()]
        rows += [
            ("tracking_efficiency_charge_separated", name, data["found"], data["truth"])
            for name, data in electron_charge_counts.items()
        ]
        rows.append(("matched_track_purity", "all", purity_num, purity_den))
        if schema["reco_pdg"]:
            rows += [("pid_purity", name, pid_num[name], pid_den[name]) for name in SPECIES.values()]
        for metric, species, num, den in rows:
            value, error = binomial(num, den)
            for i in range(len(bins) - 1):
                writer.writerow([metric, species, bins[i], bins[i + 1], num[i], den[i], value[i], error[i]])

    summary = {
        "inputs": args.files,
        "schema": schema,
        "selection": {"generatorStatus": 1, "charged": True, "eta": [args.eta_min, args.eta_max], "min_association_weight": args.min_weight},
        "totals": dict(totals),
        "overall_matched_track_purity": totals["reco_matched"] / totals["reco_accepted"] if totals["reco_accepted"] else math.nan,
        "duplicate_excess_per_matched_truth": (
            totals["duplicate_excess"] / totals["matched_truth_eligible"]
            if totals["matched_truth_eligible"]
            else None
        ),
        "notes": {
            "tracking_efficiency": "eligible status-1 charged MC particles with >=1 associated accepted reconstructed charged particle",
            "matched_track_purity": "accepted reconstructed charged particles associated to an eligible status-1 charged MC particle",
            "pid_purity": "among accepted reconstructed particles assigned a PDG species, fraction matched to the same true species",
            "duplicate_excess": "number of reconstructed particles beyond the first associated to each eligible truth particle",
        },
    }
    with (args.output / "summary.json").open("w") as stream:
        json.dump(summary, stream, indent=2, allow_nan=True)
    with (args.output / "pt_response.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["true_pt_low", "true_pt_high", "reco_pt_low", "reco_pt_high", "count"])
        for i in range(len(bins) - 1):
            for j in range(len(bins) - 1):
                writer.writerow([bins[i], bins[i + 1], bins[j], bins[j + 1], pt_response[i, j]])
    with (args.output / "momentum_resolution.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["residual", "species", "pt_low", "pt_high", "count", "median_bias", "sigma68", "mean", "standard_deviation"])
        for residual_name, collection in (("delta_pt_over_pt", pt_residuals), ("delta_p_over_p", momentum_residuals)):
            for species, values in collection.items():
                count, median, sigma68, mean, rms = resolution_statistics(values, args.min_resolution_entries)
                for i in range(len(bins) - 1):
                    writer.writerow([residual_name, species, bins[i], bins[i + 1], count[i], median[i], sigma68[i], mean[i], rms[i]])
        for residual_name, collection in (
            ("delta_pt_over_pt_charge_separated", electron_charge_pt_residuals),
            ("delta_p_over_p_charge_separated", electron_charge_momentum_residuals),
        ):
            for species, values in collection.items():
                count, median, sigma68, mean, rms = resolution_statistics(values, args.min_resolution_entries)
                for i in range(len(bins) - 1):
                    writer.writerow([residual_name, species, bins[i], bins[i + 1], count[i], median[i], sigma68[i], mean[i], rms[i]])
    if confusion:
        with (args.output / "pid_confusion.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["true_abs_pdg", "reco_abs_pdg", "count"])
            writer.writerows((true, reco, count) for (true, reco), count in sorted(confusion.items()))

    print(f"Analyzed {totals['events']:,} events; results are in {args.output}")


if __name__ == "__main__":
    main()

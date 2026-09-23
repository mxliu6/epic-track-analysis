# ePIC charged-track performance workflow

This is a ROOT-free first-pass analysis for EICrecon EDM4eic files. It uses the
`MCParticles`, reconstructed charged-particle, and MC/reconstruction association
collections stored in the `events` tree.

## 1. Set up on this Mac

The current machine is Apple Silicon and has `/usr/bin/python3` (Python 3.9), but
does not currently have CERN ROOT or the required Python analysis packages.

```bash
cd /Users/mliu/.codex/.chatgpt-projects/g-p-69602c8c645c8191a612ce2b2d83890b
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

ROOT is **not required** for these scripts. Install ROOT/eic-shell later if you
need EICrecon itself, C++ `RDataFrame` macros, EDM class dictionaries, or detector
software rather than just columnar analysis of the output.

## 2. Identify the file type and schema

```bash
python inspect_epic_file.py /absolute/path/to/file.root
python validate_epic_inputs.py /absolute/path/to/files/*.root
```

Expected reconstructed input has an `events` object and branches resembling:

- `MCParticles/MCParticles.PDG`, `.generatorStatus`, `.charge`, and `.momentum.*`
- `ReconstructedChargedParticles/ReconstructedChargedParticles.momentum.*`
- `_ReconstructedChargedParticleAssociations_rec/...rec.index`
- `_ReconstructedChargedParticleAssociations_sim/...sim.index`

The scripts also recognize the `ReconstructedChargedWithoutPIDParticles` naming
used in some productions. A HepMC file is generator input, and a DD4hep ROOT file
with simulated hits but no reconstructed collections must first be processed by
EICrecon before reconstructed-track efficiency/purity can be measured.

The validator catches incomplete copies by comparing the ROOT header's expected
end position with the actual file length, then checks the required collections.

## 3. Run the analysis

```bash
python analyze_track_performance.py \
  /absolute/path/to/sample1.root /absolute/path/to/sample2.root \
  --eta-min -3.5 --eta-max 3.5 \
  --pt-bins 0,0.2,0.5,1,2,3,5,7.5,10,15,20,30,50 \
  --output results/my_production
```

Outputs include efficiency and purity plots, a machine-readable JSON summary,
binned CSV metrics, a true-versus-reconstructed pT response matrix, and (when
reconstructed PDG is available) PID purity and a PID confusion table.
Matched tracks also produce `pt_resolution_vs_pt.png` and
`momentum_resolution_vs_pt.png`. Each shows the median residual (bias) and
`sigma68 = (q84 - q16)/2` versus true pT by truth species; the underlying counts,
robust widths, means, and standard deviations are written to
`momentum_resolution.csv`.

Electron and positron results are also separated using the signed truth PDG code
(`11` for electrons and `-11` for positrons). The dedicated outputs are
`electron_positron_tracking_efficiency.png`,
`electron_positron_pt_resolution_vs_pt.png`, and
`electron_positron_momentum_resolution_vs_pt.png`. The general species plots keep
the charge-combined category under the label `electron_or_positron`.
`electron_positron_pt_distributions.png` shows the stable-truth and
matched-reconstructed pT spectra separately, with the corresponding bin counts in
`electron_positron_pt_distributions.csv`.

## DIS-scattered electron identification

The truth DIS electron is identified from MC lineage: start from an incoming beam
electron (`PDG == 11`, `generatorStatus == 4`), follow the stored MC parent links,
and select a stable electron (`PDG == 11`, `generatorStatus == 1`) descended from
that beam particle. If radiation or conversions produce multiple stable electron
descendants, the highest-pT candidate is selected and the ambiguous event is
counted in `summary.json`. This avoids confusing unrelated conversion electrons
with the hard-scattered beam electron.

`dis_scattered_electron_pt.png` contains the truth and matched-reconstructed pT
spectra, `dis_scattered_electron_efficiency.png` contains its reconstruction
efficiency versus true pT,
and `dis_scattered_electron_pt.csv` contains the bin counts.

`dis_scattered_electron_eta.png` contains the truth and matched-reconstructed eta
spectra. `dis_scattered_electron_eta_vs_pt.png` contains side-by-side accepted
truth and matched-reconstructed 2-D eta-versus-pT distributions. Their numerical
bin contents are written to `dis_scattered_electron_eta.csv` and
`dis_scattered_electron_eta_vs_pt.csv`. Use `--eta-bins` to override the default
32 bins spanning `-4 < eta < 4`.

## Definitions

- **Truth denominator:** stable (`generatorStatus == 1`), charged MC particles in
  the requested eta range, split by absolute PDG: 11, 13, 211, 321, and 2212.
- **Tracking efficiency:** fraction of the truth denominator with at least one
  associated accepted reconstructed charged particle, binned in truth pT.
- **Matched-track purity:** fraction of accepted reconstructed charged particles
  associated to an eligible truth particle, binned in reconstructed pT.
- **Duplicate excess:** reconstructed particles beyond the first associated to a
  given eligible truth particle.
- **PID purity:** kept separate from track purity; among reconstructed particles
  assigned a species, the fraction associated to the same true species.

These are operational definitions. Before publishing results, fix the physics
acceptance (eta, vertex/origin, pT threshold), association-weight rule, treatment
of secondaries, and whether the denominator is generator-stable particles or
detector-surface particles. Record the ePIC geometry, EICrecon version, generator,
beam configuration, and reconstruction configuration with every result.

## Agentic Codex loop

1. Point Codex to one representative absolute file path.
2. Codex runs the inspector and reports the actual collections and production type.
3. Codex adapts the collection mapping and selection configuration if needed.
4. Run a small event subset/single file and validate counts, association indices,
   pT/eta spectra, and several hand-inspected events.
5. Run the full file set in chunks; compare summary JSON/CSV across productions.
6. Add regression checks so changes in efficiency, fake rate, duplicates, or PID
   response beyond agreed tolerances are flagged.

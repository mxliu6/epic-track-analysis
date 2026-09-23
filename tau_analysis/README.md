# Tau and decay-daughter analysis

This workflow follows the earlier DIS analysis pattern: inspect the schema, select particles through MC ancestry, calculate observables, validate, and export plots plus numerical data. It runs in Python using Uproot, Awkward Array, NumPy and Matplotlib. CERN ROOT and the ePIC runtime are not required.

## Setup and run

From the repository directory, create an environment if needed:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Set `TAU_DATA_DIR` to your local data directory, then run exactly two files:

```bash
TAU_DATA_DIR=/path/to/BSM-pythia/Data
.venv/bin/python tau_analysis/analyze_taus.py \
  "$TAU_DATA_DIR/bsm_ceu_275x18_100K_wo_PhotonRad_ab_0000.edm4hep.root" \
  "$TAU_DATA_DIR/bsm_ceu_275x18_100K_wo_PhotonRad_ab_1000.edm4hep.root" \
  --output tau_analysis/results
.venv/bin/python -m unittest discover -s tau_analysis -p 'test_*.py'
```

Inputs are opened read-only; output goes only to the chosen result directory. Reads are chunked, with only the needed truth branches loaded. Selected rows are accumulated in memory, suitable for this two-file analysis. The input requires a PODIO `events` TTree, `podio_metadata` collection mappings, and the MCParticle branches listed in `FIELDS` in the script, including endpoint momenta. Other schemas may need an adapter.

The [example report](example/REPORT.md) and six figures are included for browsing on GitHub. They contain no absolute local paths. Raw simulation files, environments, caches, and generated `results/` directories are excluded from Git. Full local output includes absolute input paths for provenance; review those before sharing a generated report or summary.

GitHub Actions runs the physics-logic tests and CLI help check without private simulation files. Full data processing was validated locally; CI does not reproduce the two-file analysis.

## Agentic workflow

Give Codex this reusable task:

> Read tau_analysis/README.md. Inspect exactly the two specified simulation files without modifying them. Run tau_analysis/analyze_taus.py. Check the per-file audit and summary, investigate nonzero unresolved decay leaves or failed validation, run the focused tests, and visually inspect all output figures. Report tau and daughter counts, signed PDG identities, pT and eta distributions, and any material limitations. Do not infer reconstruction performance from MC truth. Keep the input list and selection definitions explicit, and rerun after any physics-selection change.

The stages and evidence are:

1. **Inspect:** open the highest `events` ROOT cycle and resolve the MCParticles collection ID from PODIO metadata. Record each input path, size, modification time, ROOT UUID, event count, software versions and analysis-script checksum.
2. **Identify:** select `abs(PDG)==15`, nonzero generator status, and no same-signed generator tau daughter. This removes history copies. The selection is inclusive; it does not tag only the signal/hard-process tau. All selected taus in this run have status 2. Do not assume the generic filename means one tau per event.
3. **Trace:** use PODIO daughter offsets and indices, validating collection IDs, bounds and reciprocal parent links. Remove status-0 simulation secondaries. Record immediate generator daughters separately from recursively reached status-1 descendants. Stop at status 1 even if detector daughters exist. Deduplicate descendants within each tau.
4. **Calculate:** `pT = hypot(px,py)` and `eta = asinh(pz/pT)` using production-vertex MCParticles momenta in GeV. No transformation to `MCParticlesHeadOnFrameNoBeamFX`, acceptance cuts, event weights, or cross-section normalization. Preserve signed PDG IDs in tables; spectra overlay charge-conjugate species together. Neutrinos are included.
5. **Validate:** inspect charge and four-momentum closure, separately using production and endpoint tau momenta. Report unresolved generator leaves and nonfinite eta explicitly. Verify histogram counts and overflow bookkeeping. Reject broken links rather than guessing ancestry.
6. **Review:** read `results/REPORT.md`, inspect the six figures, and use the CSVs to reproduce them. The result is MC truth PID and kinematics, not detector PID or reconstructed tau efficiency.

## Outputs

- `REPORT.md`: readable findings and figures.
- `tau_distributions.png`, `direct_distributions.png`, `stable_distributions.png`: PID, pT and eta.
- `*_eta_vs_pt.png`: 2-D eta versus pT for all species at each level.
- `particles.csv`: selected particle rows with file ID, local event entry, tau index, particle index, PDG, status, charge, pT and eta.
- `decays.csv`: tau decay composition and closure diagnostics.
- `pid_counts.csv`, `decay_channels.csv`: signed identities and observed channel fractions.
- `histograms.csv`: signed-species and inclusive 1-D bin counts with square-root count errors.
- `eta_vs_pt.csv`, `histogram_accounting.csv`: 2-D counts and 1-D accounting, including nonfinite values and flow bins. The last histogram bin includes its upper edge.
- `file_audit.csv`, `summary.json`: per-file counts, provenance, definitions and validation.

Square-root count errors are descriptive Poisson errors; particles from one decay are correlated. These outputs do not provide event-level covariance or a branching-fraction measurement. Per-file event entry is the provenance identifier, not a globally unique physics event number. Duplicate input paths are rejected; cross-file physics-event duplication is not established by this analysis.

## References

- [Uproot documentation](https://uproot.readthedocs.io/en/stable/basic.html): ROOT-file reading in Python.
- [EDM4hep schema](https://github.com/key4hep/EDM4hep/blob/main/edm4hep.yaml): MCParticle momentum, status and parent/daughter relationships.

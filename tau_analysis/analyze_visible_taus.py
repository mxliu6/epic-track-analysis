#!/usr/bin/env python3
"""Visible tau observables from truth-associated reconstructed detector candidates."""
import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from analyze_taus import FIELDS, csvwrite, relations, stable_descendants
import awkward as ak
import numpy as np
import uproot
import matplotlib.pyplot as plt
from daughter_multiplicity import make_multiplicity

NEUTRINOS = {12, 14, 16}
RECO = 'ReconstructedParticles'
ASSOC = 'ReconstructedParticleAssociations'


def observables(fours, tau_p):
    """Four-vectors use (E, px, py, pz), GeV; undefined mass/angle -> NaN."""
    if not fours:
        return dict(mass_GeV=math.nan, mass2_GeV2=math.nan, scalar_pt_GeV=0.,
                    vector_pt_GeV=0., angle_deg=math.nan, energy_GeV=0.)
    a = np.asarray(fours, dtype=float)
    total = a.sum(axis=0)
    m2 = float(total[0]**2 - total[1:] @ total[1:])
    tolerance = 1e-10 * max(1., total[0]**2, float(total[1:] @ total[1:]))
    mass = math.sqrt(max(0., m2)) if m2 >= -tolerance else math.nan
    tau_p = np.asarray(tau_p)
    norm = np.linalg.norm(total[1:]) * np.linalg.norm(tau_p)
    angle = math.degrees(math.atan2(np.linalg.norm(np.cross(tau_p, total[1:])),
                                   float(tau_p @ total[1:]))) if norm > 0 else math.nan
    return dict(mass_GeV=mass, mass2_GeV2=m2,
                scalar_pt_GeV=float(np.hypot(a[:, 1], a[:, 2]).sum()),
                vector_pt_GeV=float(np.hypot(total[1], total[2])),
                angle_deg=angle, energy_GeV=float(total[0]))


def choose_matches(rec, sim, weights, charges):
    """Best truth per reco, then one reco per truth; prefer a track over a shower."""
    by_reco = defaultdict(list)
    for r, s, w in zip(rec, sim, weights):
        if math.isfinite(w) and w > 0:
            by_reco[r].append((w, s))
    by_truth = defaultdict(list)
    for r, matches in by_reco.items():
        w, s = sorted(matches, key=lambda x: (-x[0], x[1]))[0]
        by_truth[s].append((r, w))
    selected = {}
    rejected = []
    for s, candidates in by_truth.items():
        ordered = sorted(candidates, key=lambda x: (charges[x[0]] == 0, -x[1], x[0]))
        selected[s] = ordered[0]
        rejected.extend((s, r, w) for r, w in ordered[1:])
    return selected, rejected


def reco_four(momentum, energy, mass, charge):
    p = np.asarray(momentum, dtype=float)
    if not np.all(np.isfinite(p)) or not math.isfinite(mass) or mass < 0:
        return None
    if charge != 0:
        # Track momentum and the reconstruction's mass hypothesis; never truth mass.
        return np.r_[math.sqrt(float(p @ p) + mass**2), p]
    if not math.isfinite(energy) or energy < mass or np.linalg.norm(p) == 0:
        return None
    # Neutral calorimeter energy and stored direction, with stored mass hypothesis.
    return np.r_[energy, math.sqrt(max(0., energy**2 - mass**2)) * p / np.linalg.norm(p)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('files', nargs='+', type=Path)
    ap.add_argument('--output', type=Path, default=Path('tau_analysis/results/visible'))
    ap.add_argument('--step-size', type=int, default=100)
    args = ap.parse_args()
    if args.step_size < 1: ap.error('--step-size must be positive')
    files = [p.resolve() for p in args.files]
    if len(files) != len(set(files)): ap.error('Duplicate input files')
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    rows, candidates, coverage, audit, manifests = [], [], [], [], []
    totals = Counter()
    for fid, path in enumerate(files):
        counts = Counter()
        stat = path.stat()
        with uproot.open(path) as f:
            t = f['events']
            meta = f['podio_metadata'].arrays(filter_name='events___CollectionTypeInfo.*', how=dict)
            ids = dict(zip(ak.to_list(meta['events___CollectionTypeInfo.name'][0]),
                           ak.to_list(meta['events___CollectionTypeInfo.collectionID'][0])))
            branches = {x: 'MCParticles.' + x for x in FIELDS}
            for relation in ['parents', 'daughters']:
                for field in ['index', 'collectionID']:
                    branches[relation+'.'+field] = '_MCParticles_'+relation+'.'+field
            for x in ['PDG','charge','mass','energy','momentum.x','momentum.y','momentum.z']:
                branches['reco.'+x] = RECO+'.'+x
            for side in ['rec', 'sim']:
                for field in ['index','collectionID']:
                    branches[side+'.'+field] = '_'+ASSOC+'_'+side+'.'+field
            branches['weight'] = ASSOC+'.weight'
            offset = 0
            for batch in t.iterate(list(branches.values()), step_size=args.step_size, how=dict):
                data = {k: ak.to_list(batch[v]) for k, v in branches.items()}
                for local in range(len(data['PDG'])):
                    e = {k:v[local] for k,v in data.items()}
                    event = offset+local
                    pdg, status = e['PDG'], e['generatorStatus']
                    n, nr = len(pdg), len(e['reco.charge'])
                    if any(len(e[k]) != n for k in FIELDS): raise ValueError('MC array lengths')
                    if any(len(v) != nr for k,v in e.items() if k.startswith('reco.')): raise ValueError('Reco array lengths')
                    parents = relations(e, 'parents', ids['MCParticles'])
                    daughters = relations(e, 'daughters', ids['MCParticles'])
                    for i, children in enumerate(daughters):
                        if any(i not in parents[j] for j in children) or any(i not in daughters[j] for j in parents[i]):
                            raise ValueError('Nonreciprocal ancestry')
                    daughters = [[j for j in js if status[j] != 0] for js in daughters]
                    taus = [i for i in range(n) if abs(pdg[i]) == 15 and status[i] != 0
                            and not any(pdg[j] == pdg[i] for j in daughters[i])]
                    arrays = [e[k] for k in ['rec.index','sim.index','rec.collectionID','sim.collectionID','weight']]
                    if len({len(a) for a in arrays}) != 1: raise ValueError('Association array lengths')
                    for r,s,rc,sc,w in zip(*arrays):
                        if rc != ids[RECO] or sc != ids['MCParticles'] or not 0 <= r < nr or not 0 <= s < n:
                            raise ValueError('Invalid association ObjectID')
                        if not math.isfinite(w): raise ValueError('Nonfinite association weight')
                    matched, rejected = choose_matches(e['rec.index'],e['sim.index'],e['weight'],e['reco.charge'])
                    counts['events'] += 1
                    counts['taus'] += len(taus)
                    counts['reco_candidates_all_events'] += nr
                    counts['duplicate_reco_matches_all_events'] += len(rejected)
                    counts['reco_with_multiple_truth_links'] += sum(v>1 for v in Counter(e['rec.index']).values())
                    for tau in taus:
                        stable, unresolved = stable_descendants(tau, daughters, status)
                        counts['unresolved_generator_leaves'] += len(unresolved)
                        visible = [j for j in stable if abs(pdg[j]) not in NEUTRINOS]
                        counts['neutrinos_excluded'] += len(stable)-len(visible)
                        counts['visible_truth_daughters'] += len(visible)
                        tau_p = np.array([e['momentum.'+x][tau] for x in 'xyz'])
                        tf, rf, charged, neutral, matched_truth = [], [], [], [], []
                        selected_reco = set()
                        for s in visible:
                            p = np.array([e['momentum.'+x][s] for x in 'xyz'])
                            tf.append(np.r_[math.sqrt(float(p@p)+e['mass'][s]**2),p])
                            r, w = matched.get(s, (-1, math.nan))
                            four = None if r < 0 else reco_four([e['reco.momentum.'+x][r] for x in 'xyz'],
                                    e['reco.energy'][r], e['reco.mass'][r], e['reco.charge'][r])
                            usable = four is not None
                            coverage.append(dict(file_id=fid,event_entry=event,tau_index=tau,
                                truth_index=s,truth_pdg=pdg[s],truth_charge=e['charge'][s],
                                reco_index=r,usable=int(usable),association_weight=w))
                            if not usable:
                                counts['invalid_reco_fourvectors'] += r >= 0
                                continue
                            assert r not in selected_reco
                            selected_reco.add(r)
                            rf.append(four); matched_truth.append(tf[-1])
                            is_charged = e['reco.charge'][r] != 0
                            (charged if is_charged else neutral).append(four)
                            counts['reco_charged' if is_charged else 'reco_neutral'] += 1
                            counts['charged_candidates_without_assigned_pid'] += is_charged and e['reco.PDG'][r] == 0
                            counts['calorimeter_only_charged_truth'] += not is_charged and e['charge'][s] != 0
                            candidates.append(dict(file_id=fid,event_entry=event,tau_index=tau,
                                truth_index=s,truth_pdg=pdg[s],reco_index=r,reco_pdg=e['reco.PDG'][r],
                                reco_charge=e['reco.charge'][r],reco_mass_GeV=e['reco.mass'][r],weight=w,
                                energy_GeV=four[0],px_GeV=four[1],py_GeV=four[2],pz_GeV=four[3],
                                pt_GeV=float(np.hypot(four[1],four[2]))))
                        counts['taus_with_reco'] += bool(rf)
                        counts['taus_all_visible_matched'] += len(rf)==len(visible) and bool(visible)
                        counts['visible_reco_daughters'] += len(rf)
                        # Diagnose omitted intermediate generator matches, without double counting descendants.
                        todo=list(daughters[tau]); descendants=set()
                        while todo:
                            j=todo.pop()
                            if j in descendants: continue
                            descendants.add(j)
                            if status[j] != 1: todo.extend(daughters[j])
                        intermediate = [j for j in descendants if status[j] != 1 and j in matched]
                        counts['intermediate_generator_matches_omitted'] += len(intermediate)
                        base = dict(file_id=fid,event_entry=event,tau_index=tau,tau_pdg=pdg[tau],
                            n_truth_visible=len(visible),n_reco=len(rf),n_charged=len(charged),n_neutral=len(neutral),
                            all_visible_matched=int(len(rf)==len(visible) and bool(visible)),
                            intermediate_matches_omitted=len(intermediate),tau_pt_GeV=float(np.hypot(*tau_p[:2])))
                        for level, fours in [('truth_visible',tf),('reco_visible',rf),('truth_matched_subset',matched_truth),
                                             ('reco_charged',charged),('reco_neutral',neutral)]:
                            rows.append(dict(base,level=level,n_constituents=len(fours),**observables(fours,tau_p)))
                offset += len(data['PDG'])
            assert offset == t.num_entries
            manifests.append(dict(file_id=fid,filename=path.name,path=str(path),bytes=stat.st_size,
                                  mtime_ns=stat.st_mtime_ns,root_uuid=str(f.file.uuid),events=t.num_entries))
        after=path.stat()
        assert (after.st_size,after.st_mtime_ns)==(stat.st_size,stat.st_mtime_ns)
        audit.append(dict(file_id=fid,**counts)); totals.update(counts)
        print(path.name,dict(counts),flush=True)
    if not rows: raise ValueError('No tau candidates')
    csvwrite(out/'observables.csv',rows)
    csvwrite(out/'constituents.csv',candidates,fields=['file_id','event_entry','tau_index','truth_index','truth_pdg','reco_index','reco_pdg','reco_charge','reco_mass_GeV','weight','energy_GeV','px_GeV','py_GeV','pz_GeV','pt_GeV'])
    csvwrite(out/'daughter_coverage.csv',coverage)
    csvwrite(out/'file_audit.csv',audit)
    multiplicity = make_multiplicity(candidates, rows, out)
    metrics = [('mass_GeV','Visible invariant mass [GeV]'),('scalar_pt_GeV','Scalar daughter pT sum [GeV]'),
               ('vector_pt_GeV','pT of summed daughter momentum [GeV]'),('angle_deg','Angle to truth tau [degrees]')]
    hrows, flows = [], []
    # Full ranges; no overflow is silently hidden. Reco plots require >=1 constituent.
    binning={}
    for key,label in metrics:
        vals=[r[key] for r in rows if r['n_constituents'] and math.isfinite(r[key])]
        upper=max(vals,default=1)
        binning[key]=np.linspace(0,max(1,math.ceil(upper*1.01)),61)
    labels={'truth_visible':'All visible truth daughters','reco_visible':'Measured visible candidates',
            'truth_matched_subset':'Truth of matched subset','reco_charged':'Measured charged only','reco_neutral':'Measured neutral only'}
    for filename,levels,title in [
        ('visible_observables',['truth_visible','reco_visible','truth_matched_subset'],'Visible tau: truth and reconstructed detector candidates'),
        ('charged_neutral_observables',['reco_visible','reco_charged','reco_neutral'],'Measured visible tau: charged and neutral contributions')]:
        fig,axs=plt.subplots(2,2,figsize=(12,9),layout='constrained')
        for ax,(key,label) in zip(axs.flat,metrics):
            for level in levels:
                vals=[r[key] for r in rows if r['level']==level and r['n_constituents'] and math.isfinite(r[key])]
                h,edges=np.histogram(vals,bins=binning[key]);ax.stairs(h,edges,label=f'{labels[level]} (N={len(vals)})')
            ax.set(xlabel=label,ylabel='Tau candidates / bin',yscale='log');ax.grid(alpha=.2);ax.legend(fontsize=8)
        fig.suptitle(title+'\nNeutrinos excluded; MC-assisted membership, no acceptance cuts')
        fig.savefig(out/(filename+'.png'),dpi=160);plt.close(fig)
    stats={}
    for level in labels:
        rs=[r for r in rows if r['level']==level]
        stats[level]={}
        for key,label in metrics:
            eligible=[r[key] for r in rs if r['n_constituents']]
            vals=[v for v in eligible if math.isfinite(v)]
            h,edges=np.histogram(vals,bins=binning[key])
            flow=dict(level=level,variable=key,total_taus=len(rs),empty=sum(r['n_constituents']==0 for r in rs),
                      nonfinite=len(eligible)-len(vals),in_range=int(h.sum()),
                      underflow=sum(v<edges[0] for v in vals),overflow=sum(v>edges[-1] for v in vals))
            assert flow['total_taus']==sum(flow[x] for x in ['empty','nonfinite','in_range','underflow','overflow'])
            flows.append(flow)
            for low,high,count in zip(edges[:-1],edges[1:],h):hrows.append(dict(level=level,variable=key,low=low,high=high,count=int(count)))
            stats[level][key]=dict(n=len(vals),mean=float(np.mean(vals)) if vals else None,
                median=float(np.median(vals)) if vals else None,min=min(vals) if vals else None,max=max(vals) if vals else None)
    csvwrite(out/'histograms.csv',hrows);csvwrite(out/'histogram_accounting.csv',flows)
    summary=dict(totals=dict(totals),statistics=stats,inputs=manifests,accepted_daughter_multiplicity=multiplicity,
        versions=dict(uproot=uproot.__version__,awkward=ak.__version__,numpy=np.__version__),
        scripts_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('analyze_taus.py'),Path(__file__).with_name('daughter_multiplicity.py')]},
        selection='Last-copy generator taus; status-1 descendants excluding abs(PDG)=12,14,16; MC-assisted association',
        reconstruction='ReconstructedParticles only. Charged: track p with stored reco mass. Neutral: stored E and direction with stored reco mass.',
        matching='Positive weight; best truth per reco; one reco per truth preferring charged then highest weight then lowest index. No weight probability interpretation.',
        caveats=['Not a data-only tau tagger; tau direction and daughter membership use MC truth.',
                 'Neutral candidates are not identified pi0 particles; photons are counted individually, without also adding parent pi0.',
                 'Missing/unmatched daughters are not filled using truth. Kinematic plots have no acceptance cuts; daughter multiplicity/PID plots require strict reconstructed |eta|<3.5. No event weights.',
                 'Intermediate-generator and simulation-secondary associations are omitted; final status-1 matches define the sample.',
                 'A calorimeter-only match to charged truth is retained as a measured neutral candidate under its stored mass hypothesis.',
                 'Split showers/duplicate tracks reduced to one candidate per truth; merging and misassociation can bias observables.',
                 'Visible invariant mass is not the full tau mass because neutrinos are missing.',
                 'Unassigned reconstructed PID (PDG=0) is common; stored mass hypotheses are used, never substituted with truth PID masses.'],
        validation=dict(object_ids_and_ancestry='passed',histogram_accounting='passed',
                        spacelike_sums=sum(r['n_constituents']>0 and not math.isfinite(r['mass_GeV']) for r in rows)))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    report=['# Reconstructed visible tau analysis','',f'Analyzed {totals["events"]:,} events and {totals["taus"]:,} last-copy taus in two input files.' if len(files)==2 else f'Analyzed {totals["events"]:,} events.',
        '', '**Measured momenta and energies are used; MC truth provides tau membership and the reference tau direction. This is a simulation performance study, not a data-only tau reconstruction.**','',
        '## Definitions','',
        '- Visible mass: `sqrt((sum E)^2 - |sum p|^2)`.',
        '- Scalar pT sum: `sum sqrt(px^2 + py^2)`; vector-sum pT: `sqrt((sum px)^2 + (sum py)^2)`.',
        '- Opening angle: the three-dimensional angle between the truth tau production momentum and the summed reconstructed daughter momentum, in degrees.',
        '- Charged candidates use track momentum and the stored reconstructed mass hypothesis. Neutral candidates use stored energy, direction and mass hypothesis; mass-zero candidates have p = E.',
        '- Electron-, muon- and tau-neutrinos and antineutrinos are excluded. pi0 contributes through matched final photons; its parent four-vector is not added.',
        '- No match: no reconstructed constituent. Empty systems have undefined mass/angle, zero sums in CSV, and are excluded from plotted distributions.',
        '- One combined reconstructed collection avoids counting the separate charged collection twice. Duplicate truth matches prefer tracks, then larger association weight, then lower index.',
        '', '## Counts','', '| Quantity | Count |','|---|---:|']
    report += [f'| {k} | {v} |' for k,v in totals.items()]
    report += ['', '## Reconstructed visible observables (at least one candidate)','', '| Observable | Mean | Median |','|---|---:|---:|']
    for key,label in metrics:
        st=stats['reco_visible'][key]
        report.append(f'| {label} | {st["mean"]:.5g} | {st["median"]:.5g} |' if st['mean'] is not None else f'| {label} | undefined | undefined |')
    report += ['', '![Visible observables](visible_observables.png)','', '![Charged and neutral](charged_neutral_observables.png)','', '## Limitations','']
    report += ['- '+x for x in summary['caveats']]
    report += ['', '## Accepted daughter multiplicity and PID', '', 'The following plots require −3.5 < reconstructed daughter eta < 3.5. Neutrinos remain excluded. Count is per tau, not per event; zero-daughter taus are retained. Photons from pi0 decays count individually. The earlier kinematic plots keep their original selection.', '', '| Accepted daughters | Taus |', '|---:|---:|']
    report += [f'| {n} | {v} |' for n,v in multiplicity['multiplicities'].items()]
    report += ['', '![Daughter counts](daughter_multiplicity.png)', '', '![Reconstructed PID by count](daughter_pid_by_multiplicity_reconstructed.png)', '', '![Truth identity by count](daughter_pid_by_multiplicity_truth.png)', '', 'PID plots count particles in each exact multiplicity category. PDG 0 means unassigned, not a neutrino. Truth identity is shown for the same accepted reconstructed candidates; it is not detector PID performance. Each PID panel sums to multiplicity times the number of taus in its category.']
    report += ['', '## Validation','', 'Association collection IDs and bounds, reciprocal MC ancestry, and histogram accounting passed. Undefined mass/angle cases are explicitly counted. Focused unit tests cover four-vector sums and duplicate matching.', '', '## Inputs','']
    report += ['- `'+str(p)+'`' for p in files]
    (out/'REPORT.md').write_text('\n'.join(report)+'\n')
    print(json.dumps(totals,indent=2))


if __name__=='__main__': main()

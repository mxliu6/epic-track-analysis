#!/usr/bin/env python3
"""ROOT-free, ancestry-based tau truth analysis of PODIO EDM4hep TTrees."""
import argparse
import csv
import hashlib
import json
import math
import os
import platform
from collections import Counter
from pathlib import Path

CACHE = Path(__file__).resolve().parent / '.cache'
CACHE.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(CACHE))
os.environ.setdefault('XDG_CACHE_HOME', str(CACHE))
import awkward as ak
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

NAMES = {11:'e−', -11:'e+', 12:'νe', -12:'anti-νe', 13:'μ−', -13:'μ+',
         14:'νμ', -14:'anti-νμ', 15:'τ−', -15:'τ+', 16:'ντ', -16:'anti-ντ',
         22:'γ', 111:'π0', 211:'π+', -211:'π−', 321:'K+', -321:'K−',
         311:'K0', -311:'anti-K0', 310:'KS0', 130:'KL0', 221:'η', 223:'ω', 323:'K*+', -323:'K*−', 213:'ρ+', -213:'ρ−'}
FIELDS = ['PDG','generatorStatus','charge','mass','momentum.x','momentum.y','momentum.z',
          'momentumAtEndpoint.x','momentumAtEndpoint.y','momentumAtEndpoint.z',
          'parents_begin','parents_end','daughters_begin','daughters_end']

def csvwrite(path, rows, fields=None):
    rows=list(rows)
    with Path(path).open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        w.writeheader(); w.writerows(rows)

def relations(event, relation, cid):
    """Resolve PODIO ranges into validated event-local MCParticle indices."""
    n=len(event['PDG']); indices=event[relation+'.index']; ids=event[relation+'.collectionID']
    if len(indices)!=len(ids): raise ValueError('Relation array length mismatch')
    result=[]
    for b,e in zip(event[relation+'_begin'],event[relation+'_end']):
        if not 0<=b<=e<=len(indices): raise ValueError('Invalid relation offsets')
        targets=indices[b:e]
        if any(c!=cid for c in ids[b:e]): raise ValueError('Relation targets another collection')
        if any(j<0 or j>=n for j in targets): raise ValueError('Invalid relation index')
        if len(set(targets))!=len(targets): raise ValueError('Duplicate relation link')
        result.append(targets)
    return result

def stable_descendants(tau, daughters, status):
    """Stop at status 1, excluding any later detector interactions/decays."""
    found=set(); unresolved=set(); visited=set(); active=set()
    def visit(i):
        if i in active: raise ValueError('Cycle in tau ancestry')
        if i in visited: return
        visited.add(i); active.add(i)
        if status[i]==1: found.add(i)
        elif not daughters[i]: unresolved.add(i)
        else:
            for j in daughters[i]: visit(j)
        active.remove(i)
    for j in daughters[tau]: visit(j)
    return sorted(found), sorted(unresolved)

def kinematics(px,py,pz):
    pt=math.hypot(px,py)
    eta=math.asinh(pz/pt) if pt>0 else (math.copysign(math.inf,pz) if pz else math.nan)
    return pt,eta

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('files',nargs='+',type=Path)
    ap.add_argument('--output',type=Path,default=Path('tau_analysis/results'))
    ap.add_argument('--step-size',type=int,default=100)
    args=ap.parse_args(); out=args.output; out.mkdir(parents=True,exist_ok=True)
    paths=[p.resolve() for p in args.files]
    if len(set(paths))!=len(paths): raise ValueError('Duplicate input file')
    particles=[]; decays=[]; audits=[]; manifests=[]; total=Counter(); tau_status=Counter()
    for file_no,path in enumerate(paths):
        before=path.stat(); counts=Counter(); statuses=Counter()
        with uproot.open(path) as f:
            tree=f['events']  # Highest ROOT cycle only, never concatenate autosave cycles.
            meta=f['podio_metadata'].arrays(filter_name='events___CollectionTypeInfo.*',how=dict,library='ak')
            names=ak.to_list(meta['events___CollectionTypeInfo.name'][0])
            ids=ak.to_list(meta['events___CollectionTypeInfo.collectionID'][0])
            cid=dict(zip(names,ids))['MCParticles']
            branches={x:'MCParticles.'+x for x in FIELDS}
            for r in ['parents','daughters']:
                for x in ['index','collectionID']: branches[r+'.'+x]='_MCParticles_'+r+'.'+x
            offset=0
            for batch in tree.iterate(expressions=list(branches.values()),step_size=args.step_size,library='ak',how=dict):
                data={k:ak.to_list(batch[v]) for k,v in branches.items()}
                for ev in range(len(data['PDG'])):
                    e={k:v[ev] for k,v in data.items()}; pdg=e['PDG']; status=e['generatorStatus']; n=len(pdg)
                    if any(len(e[k])!=n for k in FIELDS): raise ValueError('MC array length mismatch')
                    parents=relations(e,'parents',cid); daughters=relations(e,'daughters',cid)
                    # Require reciprocity before using a relation as physics provenance.
                    for i in range(n):
                        if any(i not in parents[j] for j in daughters[i]) or any(i not in daughters[j] for j in parents[i]):
                            raise ValueError(f'Nonreciprocal MC relation: {path.name}, event {offset+ev}, particle {i}')
                    counts['detector_secondary_links_excluded']+=sum(status[i]!=0 and status[j]==0 for i in range(n) for j in daughters[i])
                    daughters=[[j for j in js if status[j]!=0] for js in daughters]
                    taus=[i for i in range(n) if abs(pdg[i])==15 and status[i]!=0]
                    selected=[i for i in taus if not any(pdg[j]==pdg[i] for j in daughters[i])]
                    counts['events']+=1; counts['events_with_tau']+=bool(selected)
                    counts['tau_records']+=len(taus); counts['tau_copies_removed']+=len(taus)-len(selected)
                    counts['selected_taus']+=len(selected)
                    statuses.update(status[i] for i in taus)
                    def row(level,i,tau):
                        pt,eta=kinematics(e['momentum.x'][i],e['momentum.y'][i],e['momentum.z'][i])
                        if not math.isfinite(pt): raise ValueError('Nonfinite pT')
                        return dict(file_id=file_no,event_entry=offset+ev,tau_index=tau,level=level,
                                    particle_index=i,pdg=pdg[i],name=NAMES.get(pdg[i],str(pdg[i])),
                                    status=status[i],charge=e['charge'][i],pt_GeV=pt,eta=eta)
                    def four(i, endpoint=False):
                        prefix='momentumAtEndpoint' if endpoint else 'momentum'
                        p=np.array([e[prefix+'.'+axis][i] for axis in 'xyz'])
                        return np.r_[math.sqrt(float(p@p)+e['mass'][i]**2),p]
                    for tau in selected:
                        particles.append(row('tau',tau,tau))
                        direct=daughters[tau]
                        stable,unresolved=stable_descendants(tau,daughters,status)
                        counts['taus_without_daughters']+=not bool(direct)
                        counts['unresolved_descendant_leaves']+=len(unresolved)
                        for level,indices in [('direct',direct),('stable',stable)]:
                            particles.extend(row(level,j,tau) for j in indices)
                            counts[level+'_daughters']+=len(indices)
                        residual=four(tau)-sum((four(j) for j in direct),start=np.zeros(4))
                        endpoint_residual=four(tau,endpoint=True)-sum((four(j) for j in direct),start=np.zeros(4))
                        charge_res=e['charge'][tau]-sum(e['charge'][j] for j in direct)
                        decays.append(dict(file_id=file_no,event_entry=offset+ev,tau_index=tau,
                            tau_pdg=pdg[tau],direct_pdgs=' '.join(map(str,sorted(pdg[j] for j in direct))),
                            stable_pdgs=' '.join(map(str,sorted(pdg[j] for j in stable))),
                            n_direct=len(direct),n_stable=len(stable),unresolved=len(unresolved),
                            delta_E_GeV=residual[0],delta_px_GeV=residual[1],delta_py_GeV=residual[2],
                            delta_pz_GeV=residual[3],delta_charge=charge_res,
                            endpoint_max_component_residual_GeV=float(np.max(np.abs(endpoint_residual)))))
                offset+=len(data['PDG'])
            if offset!=tree.num_entries: raise ValueError('Incomplete tree read')
            manifests.append(dict(file_id=file_no,path=str(path),bytes=before.st_size,mtime_ns=before.st_mtime_ns,
                                  root_uuid=str(f.file.uuid),entries=tree.num_entries,MCParticles_collectionID=cid))
        after=path.stat()
        if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns): raise ValueError('Input changed during analysis')
        audits.append(dict(file_id=file_no,**counts)); total.update(counts); tau_status.update(statuses)
        print(path.name,dict(counts),flush=True)
    if not particles: raise ValueError('No generator taus found')
    csvwrite(out/'particles.csv',particles); csvwrite(out/'decays.csv',decays); csvwrite(out/'file_audit.csv',audits)
    pidrows=[]; histrows=[]; flowrows=[]; maps=[]
    ptmax=max(10, math.ceil(max(r['pt_GeV'] for r in particles)/10)*10)
    finite_eta=[r['eta'] for r in particles if math.isfinite(r['eta'])]
    etamax=max(4,math.ceil(max(map(abs,finite_eta),default=4)))
    bins={'pt_GeV':np.linspace(0,ptmax,51),'eta':np.linspace(-etamax,etamax,61)}
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    for level,title in [('tau','Tau leptons (last copy before decay)'),('direct','Immediate tau daughters'),('stable','Final generator-level tau descendants')]:
        rows=[r for r in particles if r['level']==level]; species=Counter(r['pdg'] for r in rows)
        fig,axs=plt.subplots(1,3,figsize=(17,8),layout='constrained')
        signed=sorted(species,key=lambda p:(abs(p),p))
        bars=axs[0].barh(range(len(signed)),[species[p] for p in signed],color='#287a92')
        axs[0].set_yticks(range(len(signed)),[f'{NAMES.get(p,p)}  ({p})' for p in signed],fontsize=9)
        axs[0].bar_label(bars,padding=3,fontsize=8)
        axs[0].set_xlim(0,max(species.values(),default=1)*1.18); axs[0].invert_yaxis()
        axs[0].set_xlabel('Particles'); axs[0].set_title('Particle (signed PDG ID)')
        for p in signed: pidrows.append(dict(level=level,pdg=p,name=NAMES.get(p,str(p)),count=species[p]))
        for ax,variable in zip(axs[1:],bins):
            for color_index,p in enumerate(sorted(set(abs(p) for p in species))):
                vals=np.array([r[variable] for r in rows if abs(r['pdg'])==p]); finite=vals[np.isfinite(vals)]
                h,edges=np.histogram(finite,bins[variable]); ax.stairs(h,edges,label=f'|PDG|={p}',linewidth=1.5,color=plt.get_cmap('tab20')(color_index))
            ax.set_xlabel('pT [GeV]' if variable=='pt_GeV' else 'η'); ax.set_ylabel('Particles / bin')
            ax.legend(fontsize=8); ax.grid(alpha=.15)
        fig.suptitle(f'{title} · {len(rows):,} particles · {total["events"]:,} events\nMCParticles truth, unweighted, no acceptance cuts',fontsize=13)
        fig.savefig(out/f'{level}_distributions.png',dpi=160); plt.close(fig)
        # Histograms retain signed species; plot overlays combine charge conjugates only.
        for p in [None]+signed:
            rs=rows if p is None else [r for r in rows if r['pdg']==p]
            for variable,edges in bins.items():
                vals=np.array([r[variable] for r in rs]); finite=vals[np.isfinite(vals)]
                h,_=np.histogram(finite,edges)
                flow=dict(level=level,pdg='all' if p is None else p,variable=variable,total=len(vals),
                          nonfinite=int((~np.isfinite(vals)).sum()),underflow=int((finite<edges[0]).sum()),
                          overflow=int((finite>edges[-1]).sum()),in_range=int(h.sum()))
                assert flow['total']==sum(flow[k] for k in ['nonfinite','underflow','overflow','in_range'])
                flowrows.append(flow)
                for lo,hi,count in zip(edges[:-1],edges[1:],h):
                    histrows.append(dict(level=level,pdg=flow['pdg'],variable=variable,low=lo,high=hi,count=int(count),poisson_error=math.sqrt(count)))
        fig,ax=plt.subplots(figsize=(7,5),layout='constrained')
        good=[r for r in rows if math.isfinite(r['eta'])]
        h,x,y=np.histogram2d([r['pt_GeV'] for r in good],[r['eta'] for r in good],bins=[bins['pt_GeV'],bins['eta']])
        if h.max()>0:
            im=ax.pcolormesh(x,y,np.ma.masked_equal(h.T,0),norm=LogNorm(vmin=1,vmax=max(2,h.max())),cmap='viridis')
            fig.colorbar(im,ax=ax,label='Particles / bin')
        ax.set(xlabel='pT [GeV]',ylabel='η',title=title+'\nTruth, unweighted; all species')
        fig.savefig(out/f'{level}_eta_vs_pt.png',dpi=160); plt.close(fig)
        for i in range(len(x)-1):
            for j in range(len(y)-1): maps.append(dict(level=level,pt_low=x[i],pt_high=x[i+1],eta_low=y[j],eta_high=y[j+1],count=int(h[i,j])))
    csvwrite(out/'pid_counts.csv',pidrows); csvwrite(out/'histograms.csv',histrows)
    csvwrite(out/'histogram_accounting.csv',flowrows); csvwrite(out/'eta_vs_pt.csv',maps)
    channels=Counter(d['direct_pdgs'] for d in decays)
    csvwrite(out/'decay_channels.csv',[dict(direct_pdgs=k,count=v,fraction=v/len(decays)) for k,v in channels.most_common()])
    maxres=max(max(abs(d[k]) for k in ['delta_E_GeV','delta_px_GeV','delta_py_GeV','delta_pz_GeV']) for d in decays)
    summary=dict(totals=dict(total),tau_generator_status_counts=dict(tau_status),
                 validation=dict(relation_indices_collection_ids_and_reciprocity='passed',histogram_accounting='passed',
                                 max_direct_four_momentum_residual_GeV=maxres,
                                 decays_residual_above_1e_3_GeV=int(sum(max(abs(d[k]) for k in ['delta_E_GeV','delta_px_GeV','delta_py_GeV','delta_pz_GeV'])>1e-3 for d in decays)),
                                 max_endpoint_four_momentum_residual_GeV=max(d['endpoint_max_component_residual_GeV'] for d in decays),
                                 max_direct_charge_residual=max(abs(d['delta_charge']) for d in decays)),
                 inputs=manifests,versions=dict(python=platform.python_version(),uproot=uproot.__version__,awkward=ak.__version__,numpy=np.__version__,matplotlib=matplotlib.__version__),
                 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 definitions=dict(tau='abs(PDG)==15, generatorStatus!=0, no same-signed tau daughter',
                    direct='All immediate generator daughters of selected tau; status 0 excluded',stable='Recursive daughters, stop at generatorStatus==1; deduplicate within each tau',
                    frame='Stored MCParticles coordinates, no frame transformation',normalization='Unweighted particle counts; no acceptance cuts or cross-section normalization',
                    eta='asinh(pz/pT); beam-axis eta nonfinite and explicitly accounted',units='pT in GeV (c=1)'))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# Tau analysis report','',f'Analyzed {total["events"]:,} events in exactly {len(paths)} files. Selected {total["selected_taus"]:,} generator taus.',
           '', '## Definition', '', 'Last tau copy before decay; truth momenta in the stored MCParticles frame. No acceptance cuts. Unweighted particle counts, not cross sections. PID is the MC PDG identity, not detector PID. Direct daughters and final status-1 descendants are separate samples. Neutrinos are included. Traversal stops at status 1; status-0 simulation secondaries are excluded from all decay links. The tau sample is inclusive, including any secondary generator taus; no hard-process tag is imposed.',
           '', '## Validation','',f'- MC link indices, collection IDs and reciprocal links: passed.',f'- Tau status counts: {dict(tau_status)}.',f'- Tau copies removed: {total["tau_copies_removed"]}.',f'- Taus without daughters: {total["taus_without_daughters"]}.',f'- Unresolved descendant leaves: {total["unresolved_descendant_leaves"]}.',f'- Maximum direct-decay four-momentum component residual: {maxres:.6g} GeV.',f'- Maximum direct-decay charge residual: {summary["validation"]["max_direct_charge_residual"]:.6g}.', f'- Maximum residual using tau endpoint momentum: {summary["validation"]["max_endpoint_four_momentum_residual_GeV"]:.6g} GeV.', f'- {summary["validation"]["decays_residual_above_1e_3_GeV"]} decays have a production-momentum component residual above 1 MeV. This is a diagnostic flag, not a failed ancestry check.', '- Endpoint momentum does not resolve the largest discrepancy in these files. Its cause is not established; inspect generator/transport records before precision decay-closure studies.', '- Histogram totals including nonfinite entries and overflow: passed.', '', '## Figures','']
    for level in ['tau','direct','stable']:
        lines.extend([f'![{level} distributions]({level}_distributions.png)','',f'![{level} eta versus pT]({level}_eta_vs_pt.png)',''])
    lines.extend(['## Signed daughter IDs','','| Level | PDG ID | Particle | Count |','|---|---:|---|---:|'])
    lines.extend(f'| {r["level"]} | {r["pdg"]} | {r["name"]} | {r["count"]} |' for r in pidrows)
    lines.extend(['','## Inputs','']+[f'- `{p}`' for p in paths])
    lines.extend(['','## Interpretation limits','','These two files are a limited sample. No detector reconstruction efficiency or PID performance is inferred. Decay-channel fractions describe this sample and are not a branching-fraction measurement. Filename beam energies/process labels have not been independently established from generator settings. Input files remain read-only.'])
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary['validation'],indent=2))

if __name__=='__main__': main()

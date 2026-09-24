"""Accepted reconstructed daughter multiplicities and signed PID by multiplicity."""
import math
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
from analyze_taus import NAMES, csvwrite


def accepted_eta(px, py, pz):
    pt = math.hypot(px, py)
    eta = math.asinh(pz / pt) if pt > 0 else math.nan
    return eta, math.isfinite(eta) and -3.5 < eta < 3.5


def make_multiplicity(candidates, rows, out):
    key = lambda r: (r['file_id'], r['event_entry'], r['tau_index'])
    groups = {key(r): [] for r in rows if r['level'] == 'reco_visible'}
    accepted = []
    for c in candidates:
        eta, keep = accepted_eta(c['px_GeV'], c['py_GeV'], c['pz_GeV'])
        if keep:
            groups[key(c)].append(c)
            accepted.append(dict(c, eta=eta))
    per_tau = [dict(file_id=k[0], event_entry=k[1], tau_index=k[2],
                    n_daughters=len(v), n_charged=sum(c['reco_charge'] != 0 for c in v),
                    n_neutral=sum(c['reco_charge'] == 0 for c in v)) for k,v in groups.items()]
    counts = Counter(r['n_daughters'] for r in per_tau)
    bins = range(max(5, max(counts, default=0)) + 1)
    csvwrite(out/'accepted_daughter_counts.csv', per_tau)
    csvwrite(out/'daughter_multiplicity.csv', [dict(n_daughters=n, n_taus=counts[n]) for n in bins])
    fields = list(candidates[0]) + ['eta'] if candidates else ['file_id','event_entry','tau_index','eta']
    csvwrite(out/'accepted_constituents.csv', accepted, fields=fields)
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    bars=ax.bar(list(bins),[counts[n] for n in bins],color='#287a92')
    ax.bar_label(bars,padding=3)
    ax.set(xlabel='Number of accepted reconstructed visible daughters per tau',
           ylabel='Tau candidates',xticks=list(bins),ylim=(0,max(counts.values(),default=1)*1.15),
           title=f'Daughter multiplicity · {len(groups):,} taus\n−3.5 < reconstructed η < 3.5; neutrinos excluded')
    fig.savefig(out/'daughter_multiplicity.png',dpi=160);plt.close(fig)
    pidrows=[]
    for mode,field in [('reconstructed','reco_pdg'),('truth','truth_pdg')]:
        species=sorted({int(c[field]) for c in accepted},key=lambda p:(abs(p),p))
        fig,axes=plt.subplots(5,1,figsize=(12,16),layout='constrained')
        for n,ax in zip(range(1,6),axes):
            counter=Counter(int(c[field]) for group in groups.values() if len(group)==n for c in group)
            assert sum(counter.values())==n*counts[n]
            for p in species:
                pidrows.append(dict(pid_type=mode,n_daughters=n,n_taus=counts[n],pdg=p,
                                    name='unassigned' if p==0 else NAMES.get(p,str(p)),count=counter[p]))
            bars=ax.bar(range(len(species)),[counter[p] for p in species],color='#287a92' if mode=='reconstructed' else '#9c5c24')
            ax.bar_label(bars,padding=2,fontsize=8)
            ax.set_xticks(range(len(species)),[f'{"unassigned" if p==0 else NAMES.get(p,str(p))}\n{p}' for p in species],fontsize=9)
            ax.set(ylabel='Daughters',title=f'Exactly {n} accepted daughters · {counts[n]:,} taus · {n*counts[n]:,} daughters',
                   ylim=(0,max(counter.values(),default=1)*1.22 if counter else 1))
            if not counts[n]:ax.text(.5,.5,'No tau candidates in this category',transform=ax.transAxes,ha='center')
        fig.suptitle(('Reconstructed PID' if mode=='reconstructed' else 'MC truth identity of the same reconstructed daughters')+
                     '\nCategories use accepted reconstructed multiplicity; −3.5 < η < 3.5; no neutrinos')
        fig.savefig(out/f'daughter_pid_by_multiplicity_{mode}.png',dpi=160);plt.close(fig)
    csvwrite(out/'daughter_pid_by_multiplicity.csv',pidrows,fields=['pid_type','n_daughters','n_taus','pdg','name','count'])
    assert sum(counts.values())==len(groups)
    assert sum(n*v for n,v in counts.items())==len(accepted)
    return dict(eta_min=-3.5,eta_max=3.5,strict_boundaries=True,
                selection='Reconstructed eta of previously selected visible constituents; no tau eta or pT cut',
                n_taus=len(groups),accepted_daughters=len(accepted),rejected_daughters=len(candidates)-len(accepted),
                multiplicities={str(n):counts[n] for n in bins},pid_accounting='passed')

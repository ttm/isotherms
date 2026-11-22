#!/usr/bin/env python3
"""Produce an explicit *exploratory* historical-versus-MD plot without fitting.

This code intentionally refuses to label the 15 ps pilot as a validated prediction.
The historical series is a visually APPROXIMATED trace of Fig. 1 of
Oliveira, Taylor & Morgan, Thin Solid Films 210/211 (1992), 76-78.
It is NOT raw experimental data nor exact symbol digitization.

Usage: python scripts/compare_historical_isotherm.py --out figures/experimental_vs_md_pilot.png
Requires a locally available results/explicit_water_pilot/summary.csv.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = ROOT / 'data/oliveira_1992_fig1_approximate_trace.csv'
PILOT = ROOT / 'results/explicit_water_pilot/summary.csv'


def read_table(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline='', encoding='utf8') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f'Empty table: {path}')
    return {name: np.array([float(row[name]) for row in rows], dtype=float)
            for name in rows[0]}


def compare(exp: dict[str, np.ndarray], md: dict[str, np.ndarray]) -> dict:
    x = exp['area_a2']; y = exp['potential_mv']
    if len(x) < 5 or np.any(~np.isfinite(x)) or np.any(np.diff(x) <= 0):
        raise ValueError('Historical trace areas must be unique and increasing')
    for k in ('potential_mv', 'reading_uncertainty_mv'):
        if len(exp[k]) != len(x) or np.any(~np.isfinite(exp[k])):
            raise ValueError(f'Invalid historical {k}')
    if np.any(exp['reading_uncertainty_mv'] <= 0):
        raise ValueError('Graphical reading uncertainty must be positive')
    a = md['area_a2']; v = md['voltage_average_v']*1000
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(v)) or np.any(a < x[0]) or np.any(a > x[-1]):
        raise ValueError('Pilot areas must lie within historical trace domain')
    if len(np.unique(a)) != len(a):
        raise ValueError('Duplicate MD area')
    expected = np.interp(a, x, y)  # linear interpolation of an APPROXIMATE trace, not fitted data
    differences = v-expected
    return dict(
        experimental_provenance='manual approximate visual trace of Oliveira et al. 1992 Fig. 1; NOT raw observations',
        experimental_subphase='ultrapure water reported pH 5.6, I~2.5e-6 M, T~295 K, alpha~0.002',
        model_provenance='neutral stearic acid OpenFF/TIP3P, fixed protonation; single 5 ps + 10 ps trajectory per area',
        conditions_match=False,
        valid_model_prediction=False,
        permitted_interpretation='illustrative diagnostic of current uncalibrated, nonconverged simulation output',
        comparison=[dict(area_a2=float(area),pilot_mv=float(pred),approximately_read_experiment_mv=float(obs),
                         raw_difference_mv=float(diff))
                    for area,pred,obs,diff in zip(a,v,expected,differences)],
    )


def plot(exp: dict[str, np.ndarray], md: dict[str, np.ndarray], out: Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 13, 'axes.titlesize': 14, 'axes.labelsize': 13,
                         'xtick.labelsize': 11, 'ytick.labelsize': 11})
    result=compare(exp,md)
    x=exp['area_a2']; y=exp['potential_mv']
    area=md['area_a2']; voltage=md['voltage_average_v']*1000
    order=np.argsort(area)
    fig, (ax, ax2)=plt.subplots(1,2,figsize=(11.4,4.4),gridspec_kw={'width_ratios':[1.2,1]})
    ax.fill_between(x,y-exp['reading_uncertainty_mv'],y+exp['reading_uncertainty_mv'],
                    alpha=.14,label='Visual-reading uncertainty, not experimental error')
    ax.plot(x,y,'o-',lw=1.6,ms=3.7,label='1992 Fig. 1 approximate trace (NOT raw points)')
    ax.plot(area[order],voltage[order],'s--',lw=1.6,ms=5.2,
            label='OpenMM 15 ps PILOT (not equilibrated)')
    ax.set(xlabel=r'Area per stearic-acid molecule ($\AA^2$)',
           ylabel='Excess surface potential (mV)',xlim=(19,81),ylim=(-100,1650),
           title='Same observable, incompatible validation status')
    ax.axhline(0,color='.4',lw=.65)
    ax.legend(loc='upper right',fontsize=10,frameon=False)
    ax.grid(alpha=.15)
    ax2.plot(x,y,'o-',lw=1.7,ms=3.5)
    ax2.fill_between(x,y-exp['reading_uncertainty_mv'],y+exp['reading_uncertainty_mv'],alpha=.14)
    ax2.set(xlabel=r'Area per molecule ($\AA^2$)',ylabel='1992 graphical trace (mV)',
            xlim=(19,81),ylim=(-45,320),title='Experimental transition region (zoom)')
    ax2.axhline(0,color='.4',lw=.65)
    ax2.grid(alpha=.15)
    fig.suptitle('Historical stearic-acid ΔV–A vs nonconverged molecular forward-model output',fontsize=15)
    fig.text(.05,.012,'Graphical reconstruction only; experimental and model pH/ionization conditions differ. '
            'No agreement, onset, or mechanism is established.',fontsize=10)
    fig.subplots_adjust(left=.09,bottom=.18,right=.97,top=.81,wspace=.32)
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out,dpi=200)
    plt.close(fig)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--historical',type=Path,default=EXPERIMENT)
    p.add_argument('--pilot',type=Path,default=PILOT)
    p.add_argument('--out',type=Path,default=ROOT/'figures/experimental_vs_md_pilot.png')
    p.add_argument('--diagnostics',type=Path,default=ROOT/'results/historical_comparison_diagnostics.json')
    args=p.parse_args()
    info=plot(read_table(args.historical),read_table(args.pilot),args.out)
    args.diagnostics.parent.mkdir(parents=True,exist_ok=True)
    args.diagnostics.write_text(json.dumps(info,indent=2)+'\n')
    print('Wrote:',args.out,'and',args.diagnostics)
    print('NOT validated: approximate figure trace; MD pilot is 5+10 ps, fixed neutral acid, single run.')

if __name__=='__main__':main()
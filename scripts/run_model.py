#!/usr/bin/env python3
"""Generate the illustrative isotherm, contribution plot and CSV.

Usage: python scripts/run_model.py --out results
       python scripts/run_model.py --out results --mc --size 32 --sweeps 500
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import asdict
from pathlib import Path
import numpy as np
from isotherms.molecular import Parameters, mean_field_point, metropolis_point
from isotherms.legacy import legacy_1992_point


def main() -> None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=Path('results'))
    ap.add_argument('--mc',action='store_true',help='sample the 2D molecular lattice rather than mean field')
    ap.add_argument('--size',type=int,default=32)
    ap.add_argument('--sweeps',type=int,default=400)
    ap.add_argument('--seed',type=int,default=2026)
    ap.add_argument('--points',type=int,default=121)
    args=ap.parse_args()
    if args.points<3: ap.error('--points must be >= 3')
    args.out.mkdir(parents=True,exist_ok=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    p=Parameters()
    areas=np.linspace(80,20,args.points)
    if args.mc:
        model=[metropolis_point(float(a),p,size=args.size,burn_sweeps=args.sweeps,
                  sample_sweeps=args.sweeps,seed=args.seed+i)[0] for i,a in enumerate(areas)]
    else:
        model=[mean_field_point(float(a),p) for a in areas]
    with (args.out/'illustrative_predictions.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(asdict(model[0])))
        writer.writeheader();writer.writerows(asdict(pt) for pt in model)
    historical=[legacy_1992_point(float(a)) for a in areas]
    x=np.array([pt.area_a2 for pt in model])
    fig,ax=plt.subplots(figsize=(7.2,4.7))
    ax.plot(x,1e3*np.array([pt.potential_v for pt in model]),label='Molecular-state model')
    ax.plot(x,1e3*np.array([pt.potential_v for pt in historical]),ls=':',label='Historical 1992 model (reconstruction)')
    ax.axhline(0,lw=.7,c='gray');ax.set(xlabel=r'Area per molecule (Å$^2$)',ylabel='Surface potential (mV)')
    ax.set_xlim(80,20);ax.legend();fig.tight_layout()
    fig.savefig(args.out/'illustrative_isotherm.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.2,4.7))
    for key,label in [('head_v','Headgroup'),('tail_v','Terminal group'),('water_v','Excess hydration water'),('double_layer_v','Diffuse double layer'),('potential_v','Total')]:
        ax.plot(x,[1e3*getattr(pt,key) for pt in model],label=label,linewidth=2.2 if key=='potential_v' else 1)
    ax.axhline(0,lw=.7,c='gray');ax.set(xlabel=r'Area per molecule (Å$^2$)',ylabel='Contribution (mV)')
    ax.set_xlim(80,20);ax.legend();fig.tight_layout()
    fig.savefig(args.out/'molecular_contributions.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(figsize=(7.2,4.7))
    axes.plot(x,[pt.condensed_fraction for pt in model],label='Condensed fraction')
    axes.plot(x,[pt.water_order for pt in model],label='Excess water orientational order')
    axes.plot(x,[pt.tail_cosine for pt in model],label='Mean chain cos(tilt)')
    axes.set(xlabel=r'Area per molecule (Å$^2$)',ylabel='Molecular order parameters')
    axes.set_xlim(80,20);axes.legend();fig.tight_layout()
    fig.savefig(args.out/'molecular_order.png',dpi=170);plt.close(fig)
    print('Saved illustrative predictions in',args.out)
    print('NO experimental data was used to fit or validate these curves.')

if __name__=='__main__': main()

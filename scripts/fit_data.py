#!/usr/bin/env python3
"""Fit exactly three interpretable parameters to measured, independently sourced CSV.

Required CSV columns: area_a2,potential_mv,sigma_mv
No experimental dataset is bundled. Error bars must be physically meaningful.
"""
from __future__ import annotations
import argparse
import csv
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from isotherms.molecular import Parameters,mean_field_point


def load(path):
    with open(path,newline='') as f:
        rows=list(csv.DictReader(f))
    if len(rows)<8: raise ValueError('At least 8 experimental points required')
    try:
        a=np.array([float(r['area_a2']) for r in rows]); y=np.array([float(r['potential_mv']) for r in rows])
        s=np.array([float(r['sigma_mv']) for r in rows])
    except KeyError as e:
        raise ValueError(f'Missing required column: {e}') from e
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(y)) or not np.all(np.isfinite(s)):
        raise ValueError('Non-finite data')
    if np.min(a)<=0 or np.min(s)<=0: raise ValueError('Areas and uncertainty must be > 0')
    if np.unique(a).size<8: raise ValueError('At least 8 distinct areas required')
    return a,y,s


def fit(a,y,s):
    p0=Parameters()
    def mk(x):
        return replace(p0,critical_area_a2=float(x[0]),
            water_bias_at_80=float(x[1]),water_condensed_shift=float(x[2]))
    def fun(x):
        p=mk(x)
        return (np.array([1000*mean_field_point(ai,p).potential_v for ai in a])-y)/s
    res=least_squares(fun,[30,-.37,.09],bounds=([20,-1,-.5],[50,0,.6]),
                      x_scale='jac',max_nfev=800,ftol=1e-9,xtol=1e-9,gtol=1e-9)
    x=res.x
    # Formal local covariance, not valid if the model is misspecified or data correlated.
    covariance=None
    if np.linalg.matrix_rank(res.jac)==3 and len(a)>3:
        covariance=(np.linalg.inv(res.jac.T@res.jac)* (np.sum(res.fun**2)/(len(a)-3))).tolist()
    return mk(x),res,covariance


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv',type=Path)
    parser.add_argument('--out',type=Path,default=Path('results/fit'))
    args=parser.parse_args()
    a,y,s=load(args.csv)
    params,res,cov=fit(a,y,s)
    args.out.mkdir(parents=True,exist_ok=True)
    report=dict(success=bool(res.success),message=res.message,
      n_points=len(a),n_fit_parameters=3,chi2=float(np.sum(res.fun**2)),
      reduced_chi2=float(np.sum(res.fun**2)/(len(a)-3)),
      fitted_parameters=dict(critical_area_a2=params.critical_area_a2,
         water_bias_at_80=params.water_bias_at_80,
         water_condensed_shift=params.water_condensed_shift),
      local_covariance=cov,notes='Illustrative model; no independent validation, model uncertainty, or correlated-error treatment.')
    (args.out/'fit.json').write_text(json.dumps(report,indent=2)+'\n')
    with (args.out/'predictions.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['area_a2','observed_mv','predicted_mv','sigma_mv'])
        for aa,yy,ss in zip(a,y,s):writer.writerow([aa,yy,1000*mean_field_point(aa,params).potential_v,ss])
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()

#!/usr/bin/env python3
"""Explore, but do not validate, affine and area-low-pass changes to an MD pilot.

Inputs: the approximate 1992 Fig. 1 trace and the 5+10 ps MD pilot.
Neither output nor fitted coefficients are independent predictions.  Area filtering
is a descriptive operation on six correlated, non-equilibrated area observations,
NOT thermodynamic averaging or removal of temporal noise.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = ROOT / 'data/oliveira_1992_fig1_approximate_trace.csv'
MD = ROOT / 'results/explicit_water_pilot/summary.csv'


def _table(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline='', encoding='utf8') as f:
        data = list(csv.DictReader(f))
    if not data:
        raise ValueError(f'Empty CSV: {path}')
    return {col: np.array([float(r[col]) for r in data]) for col in data[0]}


def rms(x: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(x))))


def gaussian_area_average(area: np.ndarray, volts_mv: np.ndarray, bandwidth_a2: float) -> np.ndarray:
    """Descriptive Gaussian smoothing in area, NO scientific uncertainty correction.

    Samples are irregularly spaced. Kernel weights are normalized at each
    requested area; no extrapolation or implicit uniform-grid resampling.
    """
    if not np.isfinite(bandwidth_a2) or bandwidth_a2 <= 0:
        raise ValueError('bandwidth_a2 must be positive')
    d = (area[:, None]-area[None, :])/bandwidth_a2
    w = np.exp(-0.5*d*d)
    return (w@volts_mv)/w.sum(axis=1)


def diagnose(historical: dict[str,np.ndarray], md: dict[str,np.ndarray], bandwidth_a2: float=2.0) -> tuple[dict,dict]:
    for key in ('area_a2','potential_mv','reading_uncertainty_mv'):
        if key not in historical: raise ValueError(f'Missing historical {key}')
    for key in ('area_a2','voltage_average_v'):
        if key not in md: raise ValueError(f'Missing MD {key}')
    x, y, e = [np.asarray(historical[k],dtype=float) for k in ('area_a2','potential_mv','reading_uncertainty_mv')]
    xa, v = [np.asarray(md[k],dtype=float) for k in ('area_a2','voltage_average_v')]
    if not (len(x)>3 and len(xa)>=4 and len(x)==len(y)==len(e) and len(xa)==len(v)):
        raise ValueError('Insufficient or incompatible points')
    if not (np.all(np.isfinite(x)) and np.all(np.diff(x)>0) and np.all(e>0)
        and np.all(np.isfinite(y)) and np.all(np.isfinite(xa)) and np.all(np.isfinite(v))):
        raise ValueError('Nonfinite or unsorted data')
    if len(np.unique(xa))!=len(xa) or xa.min()<x.min() or xa.max()>x.max():
        raise ValueError('MD area outside historical domain or duplicated')
    ix=np.argsort(xa);xa=xa[ix];v=v[ix]*1000
    yi=np.interp(xa,x,y)
    visual_err=np.interp(xa,x,e)
    # Unweighted OLS avoids treating subjective visual-reading errors as rigorous
    # experimental standard deviations. This IS an empirical fit, not calibration
    # of a physical response function.
    X=np.column_stack((v,np.ones(len(v))))
    gain, offset=np.linalg.lstsq(X,yi,rcond=None)[0]
    best = gain*v + offset
    # Fit at endpoints, then test all intermediate areas. Still empirical anchoring.
    gain_end=(yi[0]-yi[-1])/(v[0]-v[-1])
    offset_end=yi[-1]-gain_end*v[-1]
    endpoints=gain_end*v+offset_end
    # A single illustrative bandwidth is chosen ex ante and never optimized.
    vs=gaussian_area_average(xa,v,bandwidth_a2)
    smoothed_calibrated=gain*vs+offset
    loo=[]
    for k in range(len(v)):
        fit=np.linalg.lstsq(np.delete(X,k,axis=0),np.delete(yi,k),rcond=None)[0]
        loo.append(float(X[k]@fit))
    loo=np.asarray(loo)
    # sign violation: curve is expected to DECREASE with increasing area in data
    wrong_intervals=[(float(xa[k]),float(xa[k+1]),float(v[k+1]-v[k]),float(yi[k+1]-yi[k]))
      for k in range(len(v)-1) if (v[k+1]-v[k])*(yi[k+1]-yi[k])<0]
    # R2 of in-sample affine fit, not held-out predictive performance.
    r2=float(1-np.sum((best-yi)**2)/np.sum((yi-np.mean(yi))**2))
    summary={
        'status':'DESCRIPTIVE CALIBRATION ONLY; nonconverged 15 ps MD and approximate experimental plot',
        'conditions_match':False,'validated_molecular_prediction':False,
        'n_md_areas':int(len(xa)),'gaussian_bandwidth_a2_illustrative':float(bandwidth_a2),
        'raw_rmse_mv':rms(v-yi),
        'fit_affine_gain':float(gain),'fit_affine_offset_mv':float(offset),
        'fit_affine_rmse_mv':rms(best-yi),'fit_affine_r2_in_sample':r2,
        'fit_affine_loocv_rmse_mv':rms(loo-yi),
        'endpoint_anchor_gain':float(gain_end),
        'endpoint_anchor_offset_mv':float(offset_end),
        'endpoint_anchor_rmse_mv_all_points':rms(endpoints-yi),
        'fixed_affine_after_gaussian_rmse_mv':rms(smoothed_calibrated-yi),
        'opposite_slope_intervals_a2':wrong_intervals,
        'notes':[
          'Offset and gain are NOT justified electrostatic corrections unless independently derived.',
          'Area-domain Gaussian smoothing cannot replace time averaging, independent replicas or phase equilibration.',
          'Historical reading uncertainty describes approximate extraction of printed dots; not laboratory measurement variance.',
          'No observation is independently predicted after the affine coefficients are fitted from all six points.'
        ]}
    data={'area_a2':xa,'md_raw_mv':v,'historical_approx_mv':yi,
       'visual_reading_uncertainty_mv':visual_err,'affine_ols_mv':best,
       'affine_anchor_mv':endpoints,'fixed_affine_after_gaussian_mv':smoothed_calibrated,
       'affine_loo_mv':loo}
    return summary,data


def plot(h:dict,m:dict,d:dict,s:dict,out:Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,2,figsize=(10.5,4.5))
    x=h['area_a2']; y=h['potential_mv']; unc=h['reading_uncertainty_mv']; a=d['area_a2']
    axs[0].fill_between(x,y-unc,y+unc,alpha=.15,label='Visual reading range, not experimental error')
    axs[0].plot(x,y,'o-',label='1992 figure (approximate)')
    axs[0].plot(a,d['md_raw_mv'],'s--',label='Unmodified 15 ps MD pilot')
    axs[0].set(xlabel=r'Area per molecule (Å$^2$)',ylabel='Excess potential (mV)',title='Unmodified values',xlim=(19,48))
    axs[0].legend(fontsize=8)
    axs[1].fill_between(x,y-unc,y+unc,alpha=.15,label='Visual reading range only')
    axs[1].plot(x,y,'o-',label='1992 figure (approximate)')
    axs[1].plot(a,d['affine_ols_mv'],'s--',label='Empirical affine fit (6 points)')
    axs[1].plot(a,d['affine_anchor_mv'],'d:',label='Anchored at 24 and 45 Å²')
    axs[1].plot(a,d['fixed_affine_after_gaussian_mv'],'^-.',alpha=.85,label='OLS coefficients + area smoothing')
    axs[1].set(xlabel=r'Area per molecule (Å$^2$)',ylabel='Excess potential (mV)',title='Diagnostically adjusted (NOT predicted)',xlim=(19,48))
    axs[1].legend(fontsize=8)
    for ax in axs:
        ax.axhline(0,lw=.6,color='gray');ax.grid(alpha=.15)
    fig.suptitle('Can scaling or smoothing reconcile a nonconverged MD pilot with the historical plot?',fontsize=12)
    fig.text(.5,.009,'Affine coefficients fitted to approximate historical points. Low-pass is area smoothing, not MD equilibration.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.045,1,.93))
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out,dpi=175);plt.close(fig)


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--historical',type=Path,default=EXPERIMENT)
    p.add_argument('--md',type=Path,default=MD)
    p.add_argument('--bandwidth-a2',type=float,default=2.)
    p.add_argument('--out',type=Path,default=ROOT/'figures/reconciliation_diagnostic.png')
    p.add_argument('--json',type=Path,default=ROOT/'results/reconciliation_diagnostic.json')
    args=p.parse_args()
    h,m=_table(args.historical),_table(args.md)
    summary,data=diagnose(h,m,args.bandwidth_a2)
    plot(h,m,data,summary,args.out)
    args.json.parent.mkdir(parents=True,exist_ok=True)
    args.json.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
#!/usr/bin/env python3
"""Compute a forward potential from atomistic coordinates, or an analytic toy benchmark.

Real input NPZ: z_nm (frames, atoms), charges_e (atoms,) OR (frames, atoms),
box_nm (3,). This is a trajectory postprocessor, not an MD simulator.

Usage:
 python scripts/atomistic_forward.py --demo --out results/atomistic
 python scripts/atomistic_forward.py --film film.npz --clean clean.npz \
     --water 5 7 --air 10 11 --out results/atomistic
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np
from isotherms.atomistic import (
    EPSILON_0, ELEMENTARY_CHARGE, NANOMETER, planar_poisson,
    surface_step, excess_surface_potential,
)


def symmetric_sheet(area_a2:float, *, n=64, q=0.012, separation_nm=.30):
    xy = (n*area_a2/100)**.5
    box = (xy,xy,14.)
    z = np.repeat([3-separation_nm/2,3+separation_nm/2,
                   9-separation_nm/2,9+separation_nm/2],n)
    charges = np.repeat([q,-q,-q,q],n)
    return planar_poisson(z, charges, box_nm=box, bins=2048), box


def save_profile(profile, path):
    arr=np.column_stack([profile.z_nm,profile.charge_density_c_m3,profile.potential_v,profile.electric_field_v_m])
    np.savetxt(path,arr,delimiter=',',comments='',header='z_nm,rho_c_m3,phi_v,field_v_m')


def demo(out:Path):
    """Demonstrate the solver, not a physical stearic-acid prediction."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    areas=np.array([20,24,28,30,35,40,50,60,80.],dtype=float)
    modeled=[]; analytic=[]
    for area in areas:
        p,box=symmetric_sheet(area)
        modeled.append(surface_step(p,bulk_water_nm=(5,7),adjacent_air_nm=(10,11)))
        S=box[0]*box[1]*NANOMETER**2
        analytic.append(64*.012*ELEMENTARY_CHARGE*.30*NANOMETER/(EPSILON_0*S))
    modeled=np.array(modeled);analytic=np.array(analytic)
    if not np.allclose(modeled,analytic,rtol=.012):
        raise RuntimeError('Toy sheet benchmark did not match analytical potential')
    with (out/'atomistic_toy_sheet.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['area_a2','periodic_poisson_v','helmholtz_analytic_v','difference_v'])
        w.writerows(zip(areas,modeled,analytic,modeled-analytic))
    p,box=symmetric_sheet(25)
    save_profile(p,out/'atomistic_toy_profile.csv')
    fig,axs=plt.subplots(1,2,figsize=(10.2,3.8))
    axs[0].plot(areas,1000*modeled,'o-',label='Atomic-sheet Poisson')
    axs[0].plot(areas,1000*analytic,'--',label='Analytical neutral sheets')
    axs[0].set(xlabel=r'Area per molecule ($\mathrm{\AA}^2$)',ylabel='Excess step (mV)',title='Code benchmark, NOT experimental data')
    axs[0].legend(fontsize=8)
    axs[1].plot(p.z_nm,p.potential_v)
    axs[1].axvspan(5,7,alpha=.13,label='bulk-water window')
    axs[1].axvspan(10,11,alpha=.13,label='air window')
    axs[1].set(xlabel='Slab coordinate z (nm)',ylabel='Electrostatic potential (V)',title='Symmetric neutral dipole sheets')
    axs[1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'atomistic_sheet_benchmark.png',dpi=180)
    plt.close(fig)
    print('Synthetic benchmark only. No water MD, stearic acid MD, or empirical fitting performed.')
    print('Max relative analytic discrepancy: %.3g' % np.max(np.abs((modeled-analytic)/analytic)))


def read_npz(file:Path,bins:int):
    with np.load(file,allow_pickle=False) as data:
        for field in ('z_nm','charges_e','box_nm'):
            if field not in data:
                raise ValueError(f'{file} missing required {field}')
        box=np.asarray(data['box_nm'],dtype=float)
        if box.shape!=(3,):
            raise ValueError(f'{file}: box_nm must be a constant 3-vector (fixed area)')
        p=planar_poisson(data['z_nm'],data['charges_e'],box_nm=tuple(box),bins=bins)
    return p,box


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--demo',action='store_true')
    p.add_argument('--film',type=Path)
    p.add_argument('--clean',type=Path)
    p.add_argument('--water',type=float,nargs=2,metavar=('START_NM','END_NM'))
    p.add_argument('--air',type=float,nargs=2,metavar=('START_NM','END_NM'))
    p.add_argument('--bins',type=int,default=1024)
    p.add_argument('--out',type=Path,default=Path('results/atomistic'))
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    if args.demo:
        if args.film or args.clean: p.error('--demo cannot be combined with --film/--clean')
        demo(args.out);return
    if not all([args.film,args.clean,args.water,args.air]):
        p.error('need --film, --clean, --water and --air for real trajectory inputs')
    film,bf=read_npz(args.film,args.bins)
    clean,bc=read_npz(args.clean,args.bins)
    if not np.allclose(bf,bc,rtol=0,atol=1e-10):
        p.error('film and clean-water trajectories require identical box dimensions; adjust reference simulation')
    excess=excess_surface_potential(film,clean,bulk_water_nm=tuple(args.water),adjacent_air_nm=tuple(args.air))
    save_profile(film,args.out/'film_profile.csv')
    save_profile(clean,args.out/'clean_profile.csv')
    with (args.out/'result.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['delta_v','film_step_v','clean_step_v','Lx_nm','Ly_nm','Lz_nm'])
        writer.writerow([excess,surface_step(film,bulk_water_nm=tuple(args.water),adjacent_air_nm=tuple(args.air)),
                         surface_step(clean,bulk_water_nm=tuple(args.water),adjacent_air_nm=tuple(args.air)),*bf])
    print(f'Excess potential: {excess*1000:+.4f} mV')
    print('Check plateau flatness, electrostatic boundary correction, finite size and replicate/block errors before interpretation.')

if __name__=='__main__':main()

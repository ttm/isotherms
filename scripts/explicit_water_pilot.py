#!/usr/bin/env python3
"""Neutral stearic-acid / explicit-TIP3P Langmuir slab research pilot.

Uses OpenFF 2.2.1 (small molecule) with AM1-BCC charges, OpenMM, and TIP3P.
Produces independent, area-constrained NVT trajectories for 2 symmetric monolayers
and a matched bare-water control at each area. Does NOT fit onset or voltage.

Install: conda/micromamba -c conda-forge openmm openmmforcefields openff-toolkit rdkit numpy scipy matplotlib
Example: python scripts/explicit_water_pilot.py --areas 24 28 30 32 36 45 --equil-ps 5 --prod-ps 10 --out results/md_pilot
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import sys
from pathlib import Path
from time import perf_counter

import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms
except ImportError as exc:
    raise RuntimeError("RDKit is required; install the documented conda environment") from exc

# Geometry deliberately independent of OpenMM so local / CI tests can inspect it.
SMILES = "CCCCCCCCCCCCCCCCCC(=O)O"  # C18H36O2, neutral acid
WATER_OH_NM = 0.09572
WATER_HOH_DEG = 104.52


def stearic_conformer(seed: int = 41):
    mol = Chem.AddHs(Chem.MolFromSmiles(SMILES))
    if Chem.rdMolDescriptors.CalcMolFormula(mol) != "C18H36O2":
        raise AssertionError("Incorrect stearic-acid molecule")
    status = AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=False)
    if status != 0:
        raise RuntimeError(f"RDKit conformer generation failed: {status}")
    # All-trans backbone, defined by the chemical graph, not a selected area.
    chain = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C"]
    for i in range(len(chain)-3):
        rdMolTransforms.SetDihedralDeg(mol.GetConformer(),*chain[i:i+4],180.0)
    ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
    ff.Minimize(maxIts=500)
    x = np.array(mol.GetConformer().GetPositions()) * 0.1  # A -> nm
    head_carbon = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C' and
                       sum(n.GetSymbol()=='O' for n in a.GetNeighbors()) == 2)
    tail_carbon = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol()=='C' and
                       sum(n.GetSymbol()=='C' for n in a.GetNeighbors()) == 1 and a.GetIdx() != head_carbon)
    head = np.mean([x[a.GetIdx()] for a in mol.GetAtoms() if a.GetSymbol()=='O'],axis=0)
    axis = x[tail_carbon] - x[head_carbon]
    axis /= np.linalg.norm(axis)
    # Rotate axis to +z. Handles numerical degeneracy near antiparallel vectors.
    up = np.array([0.,0.,1.])
    c = np.dot(axis,up)
    if c < -0.999999:
        rot = np.diag([1., -1., -1.])
    elif c > 0.999999:
        rot = np.eye(3)
    else:
        v = np.cross(axis,up)
        skew = np.array([[0.,-v[2],v[1]],[v[2],0.,-v[0]],[-v[1],v[0],0.]])
        rot = np.eye(3) + skew + skew@skew/(1.+c)
    x = (x-head) @ rot.T
    return mol,x


def _spin_z(theta):
    return np.array([[np.cos(theta),-np.sin(theta),0.],
                     [np.sin(theta), np.cos(theta),0.],[0.,0.,1.]])


def make_coordinates(area_a2: float, *, side=4, seed=41, water_thickness_nm=4.0, box_z_nm=16.0,
                     water_density_nm3=33.0):
    """Return atomic positions in nm for 2x(side**2) acids + explicit H2O slab.

    The model has *no empirical onset parameter*; A fixes only periodic xy size.
    Neutral COOH is an acidic-subphase limiting case; not constant-pH dynamics.
    Water lattice is an initial condition, equilibrated by OpenMM before analysis.
    """
    if side < 2 or area_a2 <= 15 or water_thickness_nm < 2 or box_z_nm < water_thickness_nm + 7:
        raise ValueError('Unphysical pilot geometry')
    mol,base=stearic_conformer(seed)
    rng=np.random.default_rng(seed)
    Lxy=side*np.sqrt(area_a2*0.01)
    box=np.array([Lxy,Lxy,box_z_nm])
    center_z=box_z_nm/2
    low_water=center_z-water_thickness_nm/2
    high_water=center_z+water_thickness_nm/2
    positions=[];elements=[];bonds=[];residues=[]
    # per-interface outward-pointing acids: lower -z, upper +z.
    step=Lxy/side
    for leaflet in (0,1):
        for i in range(side):
            for j in range(side):
                phi=rng.uniform(0,2*np.pi)
                # initial tilts are randomized but NOT selected to reproduce a voltage transition
                tilt=rng.uniform(0,0.12)
                ry=np.array([[np.cos(tilt),0,np.sin(tilt)],[0,1,0],[-np.sin(tilt),0,np.cos(tilt)]])
                xyz=base @ (ry @ _spin_z(phi)).T
                if leaflet==0: xyz[:,2]*=-1
                xyz += np.array([(i+.5)*step,(j+.5)*step,low_water-0.07 if leaflet==0 else high_water+0.07])
                start=len(positions)
                positions.extend(xyz)
                elements.extend([a.GetSymbol() for a in mol.GetAtoms()])
                bonds.extend((start+b.GetBeginAtomIdx(), start+b.GetEndAtomIdx()) for b in mol.GetBonds())
                residues.append(('SA',start,len(positions)))
    n_acid_atoms=len(positions)
    # Near-density-matched simple cubic O lattice for the initial water slab.
    d=water_density_nm3**(-1/3)
    nx=max(2,int(round(Lxy/d))); ny=nx
    nz=max(2,int(round(water_thickness_nm/d)))
    xs=(np.arange(nx)+.5)*Lxy/nx
    ys=(np.arange(ny)+.5)*Lxy/ny
    zs=low_water+(np.arange(nz)+.5)*water_thickness_nm/nz
    # TIP3P water: O-H distances and angle, random rigid rotations.
    t=math.radians(WATER_HOH_DEG)
    water=np.array([[0,0,0],[WATER_OH_NM,0,0],[WATER_OH_NM*math.cos(t),WATER_OH_NM*math.sin(t),0]])
    for x in xs:
        for y in ys:
            for z in zs:
                u=rng.normal(size=(3,3));Q,_=np.linalg.qr(u)
                if np.linalg.det(Q)<0: Q[:,0]*=-1
                coords=water @ Q.T + np.array([x,y,z])
                start=len(positions)
                positions.extend(coords)
                elements.extend(['O','H','H'])
                bonds.extend([(start,start+1),(start,start+2)])
                residues.append(('HOH',start,start+3))
    return {
        'mol':mol,'positions_nm':np.array(positions), 'elements':elements, 'bonds':bonds,
        'residues':residues, 'box_nm':box,'n_acid_per_leaflet':side**2,
        'n_acid_atoms':n_acid_atoms,'n_water':nx*ny*nz,
        'water_o_indices':np.arange(n_acid_atoms,len(positions),3,dtype=int),
        'surface_z_nm':(low_water,high_water),
    }


def make_openmm_topology(geom, include_acid=True):
    from openmm.app import Topology, element
    from openmm import Vec3
    from openmm import unit
    top=Topology(); chain=top.addChain()
    lookup={}; positions=[]; global_ids=[]
    for res_name,start,end in geom['residues']:
        if not include_acid and res_name=='SA': continue
        res=top.addResidue(res_name,chain)
        for j in range(start,end):
            atom_name=f"{geom['elements'][j]}{j-start+1}" if res_name=='SA' else ('O','H1','H2')[j-start]
            lookup[j]=top.addAtom(atom_name,element.getBySymbol(geom['elements'][j]),res)
            global_ids.append(j)
            positions.append(geom['positions_nm'][j])
    for i,j in geom['bonds']:
        if i in lookup and j in lookup:top.addBond(lookup[i],lookup[j])
    Lx,Ly,Lz=geom['box_nm']
    top.setPeriodicBoxVectors((Vec3(Lx,0,0),Vec3(0,Ly,0),Vec3(0,0,Lz))*unit.nanometer)
    return top,np.asarray(positions),np.asarray(global_ids)


def simulation_run(geom, *, include_acid, equil_ps, prod_ps, interval_ps, temperature_k, seed, platform_name,
                   small_molecule_forcefield, out_npz):
    import openmm as mm
    from openmm import app, unit
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    from openff.toolkit import Molecule
    topology,positions,ids=make_openmm_topology(geom,include_acid=include_acid)
    ff=app.ForceField('amber14/tip3p.xml')
    if include_acid:
        off=Molecule.from_rdkit(geom['mol'],allow_undefined_stereo=True,hydrogens_are_explicit=True)
        generator=SMIRNOFFTemplateGenerator(molecules=[off],forcefield=small_molecule_forcefield)
        ff.registerTemplateGenerator(generator.generator)
    system=ff.createSystem(topology,nonbondedMethod=app.PME,
                 nonbondedCutoff=0.9*unit.nanometer,constraints=app.HBonds,rigidWater=True,
                 ewaldErrorTolerance=5e-4)
    nb=next((f for f in system.getForces() if isinstance(f,mm.NonbondedForce)),None)
    if nb is None: raise RuntimeError('No NonbondedForce, cannot get charges')
    charges=np.array([nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles())])
    if abs(np.sum(charges))>1e-5: raise RuntimeError('System not neutral; incorrectly parameterized')
    integrator=mm.LangevinMiddleIntegrator(temperature_k*unit.kelvin,1/unit.picosecond,2*unit.femtoseconds)
    integrator.setRandomNumberSeed(int(seed))
    platform=mm.Platform.getPlatformByName(platform_name)
    props={'Threads':'2'} if platform_name=='CPU' else {}
    sim=app.Simulation(topology,system,integrator,platform,props)
    sim.context.setPositions(positions*unit.nanometer)
    sim.minimizeEnergy(maxIterations=250)
    sim.context.setVelocitiesToTemperature(temperature_k*unit.kelvin,int(seed))
    steps_per_ps=500
    sim.step(max(1,int(round(equil_ps*steps_per_ps))))
    frames=[];water_o_original=set(geom['water_o_indices'])
    water_o_local=np.array([k for k,orig in enumerate(ids) if orig in water_o_original],dtype=int)
    frame_count=max(1,int(np.ceil(prod_ps/interval_ps)))
    for _ in range(frame_count):
        sim.step(max(1,int(round(interval_ps*steps_per_ps))))
        state=sim.context.getState(getPositions=True,getEnergy=True)
        xyz=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        center=float(np.mean(xyz[water_o_local,2]))
        xyz[:,2]=np.mod(xyz[:,2] +(geom['box_nm'][2]/2 - center),geom['box_nm'][2])
        frames.append(xyz[:,2].astype(np.float32))
    np.savez_compressed(out_npz,z_nm=np.array(frames),charges_e=charges,
                        box_nm=geom['box_nm'],acid_atom_mask=np.array([orig<geom['n_acid_atoms'] for orig in ids]),
                        water_oxygen_local=water_o_local,area_a2=geom['box_nm'][0]**2*100/geom['n_acid_per_leaflet'])
    return {"frames":len(frames),"n_particles":len(charges),"net_charge_e":float(np.sum(charges)),"n_water":geom['n_water']}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--areas',type=float,nargs='+',default=[24,28,30,32,36,45])
    ap.add_argument('--side',type=int,default=4)
    ap.add_argument('--equil-ps',type=float,default=5.)
    ap.add_argument('--prod-ps',type=float,default=10.)
    ap.add_argument('--interval-ps',type=float,default=1.)
    ap.add_argument('--temp-k',type=float,default=295.)
    ap.add_argument('--seed',type=int,default=2026)
    ap.add_argument('--platform',choices=['CPU','Reference','CUDA','OpenCL'],default='CPU')
    ap.add_argument('--small-molecule-forcefield',default='openff-2.2.1')
    ap.add_argument('--out',type=Path,default=Path('results/explicit_water_pilot'))
    ap.add_argument('--bins',type=int,default=512)
    args=ap.parse_args()
    if args.equil_ps<0 or args.prod_ps<=0 or args.interval_ps<=0:ap.error('bad simulation time')
    from isotherms.atomistic import planar_poisson, surface_step
    args.out.mkdir(exist_ok=True,parents=True)
    metadata={"simulation":"neutral-stea-acid-two-interface-pilot", "areas_a2":args.areas,
              "side":args.side,"equil_ps":args.equil_ps,"prod_ps":args.prod_ps,
              "interval_ps":args.interval_ps,"temp_k":args.temp_k,"seed":args.seed,
              "forcefield":args.small_molecule_forcefield,"water":"amber14/tip3p.xml",
              "integrator":"LangevinMiddle,2fs,NVT","electrostatics":"3D PME, symmetric 2-interface slab",
              "epistemic_status":"nonconverged pilot; not evidence for a transition"}
    (args.out/'protocol.json').write_text(json.dumps(metadata,indent=2)+'\n')
    summary=[];start=perf_counter()
    for idx,area in enumerate(args.areas):
        geom=make_coordinates(area,side=args.side,seed=args.seed+idx)
        name=f'A{area:g}'.replace('.','p')
        film_path=args.out/(name+'_film.npz')
        clean_path=args.out/(name+'_clean.npz')
        print('Area',area,'A2','box',geom['box_nm'],'water',geom['n_water'],flush=True)
        film=simulation_run(geom,include_acid=True,equil_ps=args.equil_ps,prod_ps=args.prod_ps,
                    interval_ps=args.interval_ps,temperature_k=args.temp_k,seed=args.seed+100*idx,
                    platform_name=args.platform,small_molecule_forcefield=args.small_molecule_forcefield,
                    out_npz=film_path)
        clean=simulation_run(geom,include_acid=False,equil_ps=args.equil_ps,prod_ps=args.prod_ps,
                    interval_ps=args.interval_ps,temperature_k=args.temp_k,seed=args.seed+100*idx+1,
                    platform_name=args.platform,small_molecule_forcefield=args.small_molecule_forcefield,
                    out_npz=clean_path)
        f=np.load(film_path); c=np.load(clean_path)
        pf=planar_poisson(f['z_nm'],f['charges_e'],box_nm=tuple(f['box_nm']),bins=args.bins)
        pc=planar_poisson(c['z_nm'],c['charges_e'],box_nm=tuple(c['box_nm']),bins=args.bins)
        # Inside water and the corresponding upper/lower air plateau (both films symmetrical).
        bw=(7.6,8.4);upper=(14.1,14.9);lower=(1.1,1.9)
        vup=(surface_step(pf,bulk_water_nm=bw,adjacent_air_nm=upper)
             -surface_step(pc,bulk_water_nm=bw,adjacent_air_nm=upper))
        vlo=(surface_step(pf,bulk_water_nm=bw,adjacent_air_nm=lower)
             -surface_step(pc,bulk_water_nm=bw,adjacent_air_nm=lower))
        # Note: two-sided agreement is a useful diagnostic but not statistical convergence.
        row={'area_a2':area,'voltage_upper_v':vup,'voltage_lower_v':vlo,
             'voltage_average_v':(vup+vlo)/2,'side_disagreement_v':abs(vup-vlo),
             'n_water':geom['n_water'],'n_acid_per_leaflet':geom['n_acid_per_leaflet'],
             'film_atoms':film['n_particles'],'clean_atoms':clean['n_particles'],
             'frames':film['frames']}
        summary.append(row)
        with (args.out/(name+'_profiles.csv')).open('w',newline='') as fout:
            w=csv.writer(fout);w.writerow(['z_nm','rho_film_c_m3','rho_clean_c_m3','phi_film_v','phi_clean_v'])
            w.writerows(zip(pf.z_nm,pf.charge_density_c_m3,pc.charge_density_c_m3,pf.potential_v,pc.potential_v))
        with (args.out/'summary.csv').open('w',newline='') as fout:
            writer=csv.DictWriter(fout,fieldnames=list(row));writer.writeheader();writer.writerows(summary)
        print(json.dumps(row),flush=True)
    metadata['elapsed_s']=round(perf_counter()-start,1)
    (args.out/'protocol.json').write_text(json.dumps(metadata,indent=2)+'\n')
    try:
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(7,4))
        ax.plot([d['area_a2'] for d in summary],[d['voltage_average_v']*1000 for d in summary], 'o-')
        ax.set(xlabel='Molecular area (Å²)',ylabel='Short-pilot excess ΔV (mV)',
               title='Exploratory explicit-water MD: not equilibrated or validated')
        fig.tight_layout();fig.savefig(args.out/'pilot_isotherm.png',dpi=160)
    except ImportError:pass
    print('PILOT ONLY: ignore apparent onset features until multi-ns, multi-size, multi-replica convergence',flush=True)


if __name__=='__main__':main()
import numpy as np
import pandas as pd
import itertools
import mdtraj as md
from mdtraj import element
import pickle
import time
import os
from argparse import ArgumentParser

def individual_rdfs(protname,path,width,proteins,residues,n_chains):
    """
    Calculate radial distribution function between the 
    centers of mass of two protein chains from a single traj. 
    width is the bin width
    of the histogram used to construct the graph.
    """
    n_chains = int(n_chains)  
    masses = []
    radiss = []
    for n in proteins.loc[protname].fasta:
        masses.append(residues.loc["p",'MW'])
        masses.append(residues.loc[n,'MW'])
        radiss.append(residues.loc["p",'sigmas']/2)
        radiss.append(residues.loc[n,'sigmas']/2)
    masses[0] += 1
    masses[-2] += 17
    # define topology that includes bead masses to calculate correct center of mass
    top = md.Topology()
    for _ in range(n_chains):
        chain = top.add_chain()
        residue = top.add_residue('C{:d}'.format(chain.index), chain, resSeq=chain.index)
        for i,resname in enumerate(proteins.loc[protname].fasta):
            # add an element with unique name to the dictionary. the letter A is prepended to avoid doubles (e.g. cysteine and carbon)
            element.Element._elements_by_symbol.pop('A'+'P', None)
            element.Element._elements_by_symbol.pop('A'+resname.upper(), None)
            el_p = element.Element.__new__(element.Element, 1, 'A'+'P', 'A'+'P', masses[2*i], radiss[2*i])
            el_n = element.Element.__new__(element.Element, 1, 'A'+resname.upper(), 'A'+resname.upper(), masses[2*i+1], radiss[2*i+1])
            atom = top.add_atom('A'+'P', element=el_p, residue=residue)
            atom = top.add_atom('A'+resname.upper(), element=el_n, residue=residue)
        for i in range(0,chain.n_atoms,2):
            if i==chain.n_atoms-2:
                top.add_bond(chain.atom(i),chain.atom(i+1))
            else:
                top.add_bond(chain.atom(i),chain.atom(i+1))
                top.add_bond(chain.atom(i),chain.atom(i+2))
    #t = md.load(path+'/{:s}.dcd'.format(protname),top=top)
    t = md.load(path+'/{:s}.dcd'.format('fix'),top=top)
    cut = 1000 # 500 ns
    t = t[cut:]
    # create trajectory and topology for centers of mass
    cmtop = md.Topology()
    cmpos = []
    for chain_ in t.top.chains:
        chain = cmtop.add_chain()
        res = cmtop.add_residue('CM', chain, resSeq=chain.index)
        cmtop.add_atom('CM', element=t.top.atom(0).element, residue=res)
        cmpos.append(md.compute_center_of_mass(
            t.atom_slice(t.top.select('chainid {:d}'.format(chain.index)))))
    cmpos = np.swapaxes(np.array(cmpos),0,1)
    cmtraj = md.Trajectory(cmpos, cmtop, t.time, t.unitcell_lengths, t.unitcell_angles)
    upper = t.unitcell_lengths[0,0] / 2 
    # calculate the rdf between the centers of mass
    for k in range(5): # divide each replica into 4 chunks
        b,e = 7800*k, 7800*(k+1)
        pairs_list = []
        for i in range(n_chains):
            for j in range(i+1,n_chains):
                pairs_list.append([i,j])
        r,rdf = md.compute_rdf(cmtraj[b:e], pairs_list, r_range=(0,upper), bin_width = width, periodic=True)
        # save results
        np.savetxt(path+'rdfs_{:d}.dat'.format(k),np.c_[r,rdf])
    
def concatenated_rdf(protname,path,width,proteins,residues,n_chains):
    """
    Caculate rdf from a single, long trajectory after concatenation 
    of several trajectories. The concatenation takes place after
    the center-of-mass trajectories have been constructed.
    The calculation will be performed for the first n_runs trajectories. 
    """
    n_chains = int(n_chains) 
    # define topology that includes bead masses to calculate correct center of mass
    masses = []
    radiss = []
    for n in proteins.loc[protname].fasta:
        masses.append(residues.loc["p",'MW'])
        masses.append(residues.loc[n,'MW'])
        radiss.append(residues.loc["p",'sigmas']/2)
        radiss.append(residues.loc[n,'sigmas']/2)
    masses[0] += 1
    masses[-2] += 17

    top = md.Topology()
    for _ in range(n_chains):
        chain = top.add_chain()
        residue = top.add_residue('C{:d}'.format(chain.index), chain, resSeq=chain.index)
        for i,resname in enumerate(proteins.loc[protname].fasta):
            # add an element with unique name to the dictionary. the letter A is prepended to avoid doubles (e.g. cysteine and carbon)
            element.Element._elements_by_symbol.pop('A'+'P', None)
            element.Element._elements_by_symbol.pop('A'+resname.upper(), None)
            el_p = element.Element.__new__(element.Element, 1, 'A'+'P', 'A'+'P', masses[2*i], radiss[2*i])
            el_n = element.Element.__new__(element.Element, 1, 'A'+resname.upper(), 'A'+resname.upper(), masses[2*i+1], radiss[2*i+1])
            atom = top.add_atom('A'+'P', element=el_p, residue=residue)
            atom = top.add_atom('A'+resname.upper(), element=el_n, residue=residue)
        for i in range(0,chain.n_atoms,2):
            if i==chain.n_atoms-2:
                top.add_bond(chain.atom(i),chain.atom(i+1))
            else:
                top.add_bond(chain.atom(i),chain.atom(i+1))
                top.add_bond(chain.atom(i),chain.atom(i+2))

    cmtrajs = []
    for run in range(1):    
        # load trajectory data 
        t = md.load(path+'/{:s}.dcd'.format("fix"),top=top)
        cut = 1000 # 100 ns
        t = t[cut:]
        # create trajectory and topology for centers of mass
        cmtop = md.Topology()
        cmpos = []
        for chain in t.top.chains:
            chain = cmtop.add_chain()
            res = cmtop.add_residue('CM', chain, resSeq=chain.index)
            cmtop.add_atom('CM', element=t.top.atom(0).element, residue=res)
            cmpos.append(md.compute_center_of_mass(
                t.atom_slice(t.top.select('chainid {:d}'.format(chain.index)))))
        cmpos = np.swapaxes(np.array(cmpos),0,1)
        cmtraj = md.Trajectory(cmpos, cmtop, t.time, t.unitcell_lengths, t.unitcell_angles)
        cmtrajs.append(cmtraj)
    # create pseudo-traj based on all individual trajectories 
    pseudo_traj = md.join(cmtrajs)
    # calculate the rdf between the centers of mass
    pairs_list = []
    for i in range(n_chains):
        for j in range(i+1,n_chains):
            pairs_list.append([i,j])
    upper = t.unitcell_lengths[0,0] / 2 
    r,rdf = md.compute_rdf(pseudo_traj, pairs_list, r_range=(0,upper), bin_width = width, periodic=True)
    # save results
    np.savetxt(path+'rdfs.dat',np.c_[r,rdf])
     

import pandas as pd
import numpy as np
import mdtraj as md
import itertools
import os
import MDAnalysis
from MDAnalysis import transformations


def calc_zpatch(z,h):
    cutoff = 20
    ct = 0.
    ct_max = 0.
    zwindow = []
    hwindow = []
    zpatch = []
    hpatch = []
    for ix, x in enumerate(h):
        if x > cutoff:
            ct += x
            zwindow.append(z[ix])
            hwindow.append(x)
        else:
            if ct > ct_max:
                ct_max = ct
                zpatch = zwindow
                hpatch = hwindow
            ct = 0.
            zwindow = []
            hwindow = []
    zpatch = np.array(zpatch)
    hpatch = np.array(hpatch)
    return zpatch, hpatch

def center_slab(path,name,temp,outname,start=None,end=None,step=1,input_pdb='top.pdb'):
    print(path)
    if not os.path.exists('traj'):
        os.system('mkdir traj')

    u = MDAnalysis.Universe(path+f'{temp}/{input_pdb}',path+f'/{temp}/{name}.dcd',in_memory=True)
    os.system(f'cp {path}/{temp}/{input_pdb} traj/{outname}_{temp}.pdb')
    n_frames = len(u.trajectory[start:end:step])
    ag = u.atoms
    n_atoms = ag.n_atoms
    print(outname, 10**-3*(u.trajectory[-1].time),"ns")
    lz = u.dimensions[2]
    edges = np.arange(0,lz+1,1)
    dz = (edges[1] - edges[0]) / 2.
    z = edges[:-1] + dz
    n_bins = len(z)
    hs = np.zeros((n_frames,n_bins))
    with MDAnalysis.Writer(f'traj/{outname}_{temp}.dcd',n_atoms) as W:
        for t,ts in enumerate(u.trajectory[start:end:step]):
            # shift max density to center
            zpos = ag.positions.T[2]
            h, e = np.histogram(zpos,bins=edges)
            zmax = z[np.argmax(h)]
            ag.translate(np.array([0,0,-zmax+0.5*lz]))
            ts = transformations.wrap(ag)(ts)
            zpos = ag.positions.T[2]
            h, e = np.histogram(zpos, bins=edges)
            zpatch, hpatch = calc_zpatch(z,h)
            zmid = np.average(zpatch,weights=hpatch)
            ag.translate(np.array([0,0,-zmid+0.5*lz]))
            ts = transformations.wrap(ag)(ts)
            zpos = ag.positions.T[2]
            h, e = np.histogram(zpos,bins=edges)
            hs[t] = h
            W.write(ag)
    if not os.path.exists('hist'):
        os.system('mkdir hist')
    np.save(f'traj/{outname}_{temp}.npy',hs,allow_pickle=False)
    return hs, z

def center_sphere(path,name,temp,outname,start=None,end=None,step=1,input_pdb='top.pdb'):
    if not os.path.exists('traj'):
        os.system('mkdir traj')

    u = MDAnalysis.Universe(path+f'{temp}/{input_pdb}',path+f'/{temp}/{name}.dcd',in_memory=True)
    os.system(f'cp {path}/{temp}/{input_pdb} traj/{outname}_{temp}.pdb')
    n_frames = len(u.trajectory[start:end:step])
    ag = u.atoms
    n_atoms = ag.n_atoms
    lz = u.dimensions[2]
    #edges = np.arange(0,lz+1,1)
    edges = np.arange(0,lz+1,1)
    dz = (edges[1] - edges[0]) / 2.
    z = edges[:-1] + dz
    n_bins = len(z)
    hs = np.zeros((n_frames,n_bins))
    with MDAnalysis.Writer(f'traj/{outname}_{temp}.dcd',n_atoms) as W:
        for t,ts in enumerate(u.trajectory[start:end:step]):
            # shift max density to center
            zpos = ag.positions.T[2]
            h, e = np.histogram(zpos,bins=edges)
            zmax = z[np.argmax(h)]
            ag.translate(np.array([0,0,-zmax+0.5*lz]))
            ts = transformations.wrap(ag)(ts)
            zpos = ag.positions.T[2]
            h, e = np.histogram(zpos, bins=edges)
            zpatch, hpatch = calc_zpatch(z,h)
            #if len(zpatch)>0:
            #    zmid = np.average(zpatch,weights=hpatch)
            #    ag.translate(np.array([0,0,-zmid+0.5*lz]))
            #    ts = transformations.wrap(ag)(ts)
            # y position
            zpos = ag.positions.T[1]
            h, e = np.histogram(zpos,bins=edges)
            zmax = z[np.argmax(h)]
            ag.translate(np.array([0,-zmax+0.5*lz,0]))
            ts = transformations.wrap(ag)(ts)
            zpos = ag.positions.T[1]
            h, e = np.histogram(zpos, bins=edges)
            zpatch, hpatch = calc_zpatch(z,h)
            #if len(zpatch)>0:
            #    zmid = np.average(zpatch,weights=hpatch)
            #    ag.translate(np.array([0,-zmid+0.5*lz,0]))
            #    ts = transformations.wrap(ag)(ts)   

            # x position
            zpos = ag.positions.T[0]
            h, e = np.histogram(zpos,bins=edges)
            zmax = z[np.argmax(h)]
            ag.translate(np.array([-zmax+0.5*lz,0,0]))
            ts = transformations.wrap(ag)(ts)
            zpos = ag.positions.T[0]
            h, e = np.histogram(zpos, bins=edges)
            zpatch, hpatch = calc_zpatch(z,h)
            #if len(zpatch)>0:
            #    zmid = np.average(zpatch,weights=hpatch)
            #    ag.translate(np.array([-zmid+0.5*lz,0,0]))
            #    ts = transformations.wrap(ag)(ts)

            W.write(ag)




#center_slab('../protein-rna/RRP','RRP',300,start=None,end=None,step=1,input_pdb='top.pdb')
#center_slab('../protein-rna/PLP','PLP',300,start=None,end=None,step=1,input_pdb='top.pdb')
#center_slab('../protein-rna_old/PLP_100_FUSRGG3_300_polyU40_50/','PLP-RRP-polyU40',300,start=None,end=None,step=1,input_pdb='top.pdb')
#center_slab('../protein-rna/FUSRGG-polyR40/FUSRGG3-120/FUSRGG3/','FUSRGG3',293,'FUSRGG3-120',start=None,end=None,step=1,input_pdb='top.pdb')
#center_slab('../protein-rna/FUSRGG-polyR40/FUSRGG3-120_polyR40-34/FUSRGG3/','FUSRGG3',293,'FUSRGG3-120_polyR40-34',start=None,end=None,step=1,input_pdb='top.pdb')
#center_slab('../protein-rna/FUSRGG-polyR40/FUSRGG3-120_polyR40-341/FUSRGG3/','FUSRGG3',293,'FUSRGG3-120_polyR40-341',start=None,end=None,step=1,input_pdb='top.pdb')


#######
# FUS
#######
cut = 20 # trj is output in every 10 ns
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-40/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-40',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-80/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-80',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-100/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-100',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-110/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-110',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-120/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-120',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-130/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-130',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-140/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-140',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-160/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-160',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-200/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-200',start=cut,end=None,step=1,input_pdb='top.pdb')

#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-40_cutoff6/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-40_cutoff6',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-80_cutoff6/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-80_cutoff6',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-100_cutoff6/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-100_cutoff6',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-110_cutoff6/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-110_cutoff6',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-120_cutoff6/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-120_cutoff6',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-130_cutoff6/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-130_cutoff6',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-140_cutoff6/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-140_cutoff6',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-160_cutoff6/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-160_cutoff6',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-polyR40_large/FUSRGG-500_polyR40-200_cutoff6/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40-200_cutoff6',start=cut,end=None,step=1,input_pdb='top.pdb')

#center_sphere('../FUSRGG-polyR40_cubic_large/FUSRGG-500/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_cubic',start=100,end=None,step=1,input_pdb='top.pdb')
#center_sphere('../FUSRGG-polyR40_cubic_large/FUSRGG-500_polyR40_62/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40_62_cubic',start=100,end=None,step=1,input_pdb='top.pdb')
#center_sphere('../FUSRGG-polyR40_cubic_large/FUSRGG-500_polyR40_750/FUSRGG3/','FUSRGG3',293,'FUSRGG3-500_polyR40_750_cubic',start=100,end=None,step=1,input_pdb='top.pdb')

#cut = 300
#center_slab('../PLP_slab/PLP-164_L150/FUSRGG3/','FUSRGG3',293,'PLP_slab',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../PLP_slab/PLP-164_L150/FUSRGG3/','FUSRGG3',297,'PLP_slab',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../PLP-polyR40_slab/PLP-164_polyU40-128/FUSRGG3/','FUSRGG3',293,'PLP-polyR40_slab',start=cut,end=None,step=1,input_pdb='top.pdb')
#center_slab('../FUSRGG-PLP-polyR40_slab/FUSRGG3/','FUSRGG3',293,'FUSRGG-PLP-polyR40_slab',start=100,end=None,step=1,input_pdb='top.pdb')
#center_sphere('../FUSRGG-PLP-polyR40_cubic/FUSRGG3/','FUSRGG3',293,'FUSRGG3-473_PLP-164_polyR40-128',start=100,end=None,step=1,input_pdb='top.pdb')



"""
MED1
output 10 ns -> 100 steps=1000 ns 
"""
step = 1 # trajectories are oputput per 10 ns
#center_slab('../MED1_slab_large/MED1-400_RNR-20/MED1/','MED1',293,'MED1-400_RNA-20',start=100,end=None,step=step,input_pdb='top.pdb') # use cutoff for patch 20
#center_slab('../MED1_slab_large/MED1-200_RNR-20/MED1/','MED1',293,'MED1-200_RNA-20',start=100,end=None,step=step,input_pdb='top.pdb')
#center_slab('../MED1_slab_large/MED1-200_RNR-40/MED1/','MED1',293,'MED1-200_RNA-40',start=100,end=None,step=step,input_pdb='top.pdb')
#center_slab('../MED1_slab_large/MED1-200_RNR-50/MED1/','MED1',293,'MED1-200_RNA-50',start=100,end=None,step=step,input_pdb='top.pdb')
#center_slab('../MED1_slab_large/MED1-200_RNR-55/MED1/','MED1',293,'MED1-200_RNA-55',start=100,end=None,step=step,input_pdb='top.pdb')
#center_slab('../MED1_slab_large/MED1-200_RNR-60/MED1/','MED1',293,'MED1-200_RNA-60',start=100,end=1000,step=step,input_pdb='top.pdb')
#center_slab('../MED1_slab_large/MED1-200_RNR-60/MED1/','MED1',293,'MED1-200_RNA-60_test',start=100,end=1500,step=20,input_pdb='top.pdb')

#center_slab('../MED1_slab_large/MED1-200_RNR-70/MED1/','MED1',293,'MED1-200_RNA-70',start=100,end=None,step=step,input_pdb='top.pdb')
#center_slab('../MED1_slab_large/MED1-200_RNR-80/MED1/','MED1',293,'MED1-200_RNA-80',start=100,end=None,step=step,input_pdb='top.pdb')
#center_slab('../MED1_slab_large/MED1-200_RNR-100/MED1/','MED1',293,'MED1-200_RNA-80',start=100,end=None,step=step,input_pdb='top.pdb')



"""
MED1+GUEST+RNA

"""
step = 1
#center_slab('../MED1_GUEST_slab_large/MED1-400_RNR-20_HP1a_200/MED1/','MED1',293,'MED1-400_RNA-20_HP1a-200',start=100,end=None,step=step,input_pdb='top.pdb')
#center_slab('../MED1_GUEST_slab_large/MED1-200_RNR-20_HP1a_100/MED1/','MED1',293,'MED1-200_RNA-20_HP1a-100',start=100,end=None,step=step,input_pdb='top.pdb')
#center_slab('../MED1_GUEST_slab_large/MED1-200_RNR-60_HP1a_100/MED1/','MED1',293,'MED1-200_RNA-60_HP1a-100',start=100,end=None,step=step,input_pdb='top.pdb')
#center_slab('../MED1_GUEST_slab_large/MED1-400_RNR-20_STP6_200/MED1/','MED1',293,'MED1-400_RNA-20_STP6-200',start=100,end=None,step=step,input_pdb='top.pdb')
#center_slab('../MED1_GUEST_slab_large/MED1-200_RNR-20_STP6_100/MED1/','MED1',293,'MED1-200_RNA-20_STP6-100',start=100,end=None,step=step,input_pdb='top.pdb')
#center_slab('../MED1_GUEST_slab_large/MED1-200_RNR-60_STP6_100/MED1/','MED1',293,'MED1-200_RNA-60_STP6-100',start=100,end=None,step=step,input_pdb='top.pdb')


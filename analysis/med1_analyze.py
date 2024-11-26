import mdtraj as md
import numpy as np
import pandas as pd
from numba import jit
import string
from scipy.ndimage import gaussian_filter1d
from mdtraj import element
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from scipy.optimize import least_squares
from scipy.stats import pearsonr, spearmanr
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler
from matplotlib.colors import LogNorm
import warnings
import itertools
warnings.filterwarnings('ignore')
import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
from statsmodels.tsa.stattools import acf
import sys
#!git clone https://github.com/fpesceKU/BLOCKING.git
sys.path.append('BLOCKING')
from main import BlockAnalysis
from PIL import Image
import matplotlib.ticker as ticker

"""
MED1 system analysis
this code calculate histgram for each component
"""




def hist_comp(name, atom_indices, n_seq,temp=293):
    trj = md.load('traj/{:s}_{:d}.dcd'.format(name,temp),top='traj/{:s}_{:d}.pdb'.format(name,temp),
                   atom_indices= atom_indices)
    z0 = trj.unitcell_lengths[0,2]
    zs = np.arange(0,z0,1) # bin size 1 nm
    pos = trj.xyz
    
    ns = []
    for t in range(trj.n_frames):
        n,_ = np.histogram(pos[t,:,2],bins=zs)
        ns.append(n)
    ns = np.array(ns)
    L = trj.unitcell_lengths[0,0]
    N = n_seq
    dLz = zs[1]-zs[0] 
    conv = 10/6.022/N/L/L/dLz*1e3 # in mM #conv = 100/6.022/N/L/L*1e3 # in mM
    xs = 0.5*(zs[1:]+zs[:-1]) 
    return xs, ns*conv


def get_csat_med1(z,h,zc_fit):
    print(h.shape)
    """Function to calculate csat and ccon from concentration profiles"""
    lz = h.shape[1]+1
    edges = np.arange(-lz/2.,lz/2.,1)/1 # bin 1 nm
    dz = (edges[1]-edges[0])/2.
    z = edges[:-1]+dz
    profile = lambda x,a,b,c,d : .5*(a+b)+.5*(b-a)*np.tanh((np.abs(x)-c)/d)
    residuals = lambda params,*args : ( args[1] - profile(args[0], *params) )
    hm = np.mean(h,axis=0)
    z1 = z[z>0]
    h1 = hm[z>0]
    z2 = z[z<0]
    h2 = hm[z<0]
    p0=[1,1,1,1]

    res1 = least_squares(residuals, x0=p0, args=[z1, h1], bounds=([0]*4,[3000]*4))
    res2 = least_squares(residuals, x0=p0, args=[z2, h2], bounds=([0]*4,[3000]*4))

    cutoffs1 = [res1.x[2]-res1.x[3],-res2.x[2]+res2.x[3]]
    cutoffs2 = [res1.x[2]+zc_fit*res1.x[3],-res2.x[2]-zc_fit*res2.x[3]]

    bool1 = np.logical_and(z<cutoffs1[0],z>cutoffs1[1])
    bool2 = np.logical_or(z>cutoffs2[0],z<cutoffs2[1])

    dilarray = np.apply_along_axis(lambda a: a[bool2].mean(), 1, h)
    denarray = np.apply_along_axis(lambda a: a[bool1].mean(), 1, h)

    dil = hm[bool2].mean()
    den = hm[bool1].mean()


    block_dil = BlockAnalysis(dilarray)
    block_den = BlockAnalysis(denarray)
    block_dil.SEM()
    block_den.SEM()
    return z,h,block_dil.av, block_dil.sem, block_den.av, block_den.sem


def get_droplet_pos(h,zc_fit):
    """
    this function return bools to indicate dense/dilute phase
    """
    lz = h.shape[1]+1
    edges = np.arange(-lz/2.,lz/2.,1)/1 # bin 1 nm
    dz = (edges[1]-edges[0])/2.
    z = edges[:-1]+dz
    profile = lambda x,a,b,c,d : .5*(a+b)+.5*(b-a)*np.tanh((np.abs(x)-c)/d)
    residuals = lambda params,*args : ( args[1] - profile(args[0], *params) )
    hm = np.mean(h,axis=0)
    z1 = z[z>0]
    h1 = hm[z>0]
    z2 = z[z<0]
    h2 = hm[z<0]
    p0=[1,1,1,1]

    res1 = least_squares(residuals, x0=p0, args=[z1, h1], bounds=([0]*4,[3000]*4))
    res2 = least_squares(residuals, x0=p0, args=[z2, h2], bounds=([0]*4,[3000]*4))

    cutoffs1 = [res1.x[2]-res1.x[3],-res2.x[2]+res2.x[3]]
    cutoffs2 = [res1.x[2]+zc_fit*res1.x[3],-res2.x[2]-zc_fit*res2.x[3]]

    bool1 = np.logical_and(z<cutoffs1[0],z>cutoffs1[1])
    bool2 = np.logical_or(z>cutoffs2[0],z<cutoffs2[1])

    return bool1, bool2


def get_conc_by_pos(h,bool1,bool2): 
    """
    this function uses bool1 and bool2 to calculate conc
    """
    dilarray = np.apply_along_axis(lambda a: a[bool2].mean(), 1, h)
    denarray = np.apply_along_axis(lambda a: a[bool1].mean(), 1, h)

    hm = np.mean(h,axis=0)
    dil = hm[bool2].mean()
    den = hm[bool1].mean()


    block_dil = BlockAnalysis(dilarray)
    block_den = BlockAnalysis(denarray)
    block_dil.SEM()
    block_den.SEM()
    return h,block_dil.av, block_dil.sem, block_den.av, block_den.sem


#========== MED1-RNA-Guest 
def main_med1_rna_hp1a(rna_num=[10,20,60]):
    guest_name = 'HP1a'
    for ind, key in enumerate(rna_num):
        [MED1_n, n_rna, guest_n] = keys[key] 
        name="MED1-{:d}_RNA-{:d}_{:s}-{:d}".format(MED1_n,n_rna,guest_name,guest_n)
        print(name)

        # MED1
        start = guest_seq[guest_name]*guest_n
        index = np.arange(start, start+MED1_seq*MED1_n)
        x1,y1 = hist_comp(name, index, MED1_seq)
        _,_,csat, csat_err, cden, cden_err = get_csat_med1(x1,y1,dil_cutoff[guest_name+"_"+str(key)])
        ## for guest
        bool_den,bool_dil = get_droplet_pos(y1,dil_cutoff[guest_name+"_"+str(key)]) 
    
        # RNA
        start = guest_seq[guest_name]*guest_n + MED1_seq*MED1_n
        index = np.arange(start, start+RNA_seq*2*n_rna)
        x2,y2 = hist_comp(name, index, RNA_seq*2)
        _,_,csat_r, csat_err_r, cden_r, cden_err_r = get_csat_med1(x2,y2,dil_cutoff[guest_name+"_"+str(key)])

        # guest
        index = np.arange(guest_seq[guest_name]*guest_n)
        x3,y3 = hist_comp(name, index, guest_seq[guest_name])
        _,csat_g, csat_err_g, cden_g, cden_err_g = get_conc_by_pos(y3,bool_den,bool_dil) 
        print("dilute",csat_g,"dense",cden_g)

        np.save("hist/PRO_MED1-{:d}_RNA-{:d}_{:s}-{:d}".format(MED1_n,n_rna,guest_name,guest_n),np.array([x1,np.mean(y1,axis=0)]))
        np.save("hist/RNA_MED1-{:d}_RNA-{:d}_{:s}-{:d}".format(MED1_n,n_rna,guest_name,guest_n),np.array([x2,np.mean(y2,axis=0)]))
        np.save("hist/GUEST_MED1-{:d}_RNA-{:d}_{:s}-{:d}".format(MED1_n,n_rna,guest_name,guest_n),np.array([x3,np.mean(y3,axis=0)]))

        np.save("fig_data/MED1-{:d}_RNA-{:d}_{:s}-{:d}.npy".format(MED1_n, n_rna,guest_name,guest_n), np.array([[csat, csat_err, cden, cden_err],
                                                                                 [csat_g, csat_err_g, cden_g, cden_err_g],
                                                                                 [csat_r, csat_err_r, cden_r, cden_err_r]]))
        print('done')


def main_med1_rna_spt6(rna_num=[10,20,60]):
    guest_name = 'STP6'
    for ind, key in enumerate(rna_num):
        [MED1_n, n_rna, guest_n] = keys[key] 
        name="MED1-{:d}_RNA-{:d}_{:s}-{:d}".format(MED1_n,n_rna,guest_name,guest_n)
        print(name)
    
        # MED1
        index = np.arange(MED1_seq*MED1_n)
        x1,y1 = hist_comp(name, index, MED1_seq)
        _,_,csat, csat_err, cden, cden_err = get_csat_med1(x1,y1,dil_cutoff[guest_name+"_"+str(key)])
        ## for guest
        bool_den,bool_dil = get_droplet_pos(y1,dil_cutoff[guest_name+"_"+str(key)])

        # RNA
        start = guest_seq[guest_name]*guest_n + MED1_seq*MED1_n
        index = np.arange(start, start+RNA_seq*2*n_rna)
        x2,y2 = hist_comp(name, index, RNA_seq*2)
        _,_,csat_r, csat_err_r, cden_r, cden_err_r = get_csat_med1(x2,y2,dil_cutoff[guest_name+"_"+str(key)])

        # guest
        start = MED1_seq*MED1_n
        index = np.arange(start, start+guest_seq[guest_name]*guest_n)
        x3,y3 = hist_comp(name, index, guest_seq[guest_name])
        _,csat_g, csat_err_g, cden_g, cden_err_g = get_conc_by_pos(y3,bool_den,bool_dil)
        print("dilute",csat_g,"dense",cden_g)

        np.save("hist/PRO_MED1-{:d}_RNA-{:d}_{:s}-{:d}".format(MED1_n,n_rna,guest_name,guest_n),np.array([x1,np.mean(y1,axis=0)]))
        np.save("hist/RNA_MED1-{:d}_RNA-{:d}_{:s}-{:d}".format(MED1_n,n_rna,guest_name,guest_n),np.array([x2,np.mean(y2,axis=0)]))
        np.save("hist/GUEST_MED1-{:d}_RNA-{:d}_{:s}-{:d}".format(MED1_n,n_rna,guest_name,guest_n),np.array([x3,np.mean(y3,axis=0)]))

        np.save("fig_data/MED1-{:d}_RNA-{:d}_{:s}-{:d}.npy".format(MED1_n, n_rna,guest_name,guest_n), np.array([[csat, csat_err, cden, cden_err],
                                                                                 [csat_g, csat_err_g, cden_g, cden_err_g],
                                                                                 [csat_r, csat_err_r, cden_r, cden_err_r]]))
        print('done')

#=========== upper section define functions

MED1_seq = 626
RNA_seq = 477
guest_seq = {'HP1a':47,'STP6':201}
dil_cutoff = {'HP1a_10':13.5,'HP1a_20':6,'HP1a_60':6,
             'STP6_10':8,'STP6_20':6,'STP6_60':6}
keys = {10:[400,20,200], 20:[200,20,100], 60:[200,60,100]}


#main_med1_rna_spt6(rna_num=[10]) #main_med1_rna_spt6(rna_num=[10,20,60])
main_med1_rna_hp1a(rna_num=[10,20,60])

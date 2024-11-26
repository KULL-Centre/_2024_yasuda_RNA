import pandas as pd
import numpy as np
import mdtraj as md
import itertools

def initProteins():
    proteins = pd.DataFrame(index=['polyA30','polyU30'], columns=['labels','eps_factor','wh','L','temp','obs','pH','ionic','expPREs','fasta','path'])
    fasta_polyA30 = """aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa""".replace('\n', '')
    fasta_polyR30 = """rrrrrrrrrrrrrrrrrrrrrrrrrrrrrr""".replace('\n', '')
    fasta_polyU30 = """uuuuuuuuuuuuuuuuuuuuuuuuuuuuuu""".replace('\n', '') 
    proteins.loc['polyA30'] = dict(labels=[16, 86, 142],eps_factor=0.2,L=50,wh=850,temp=300,obs='rate',pH=7.0,fasta=list(fasta_polyA30),ionic=0.02,path='./') 
    proteins.loc['polyR30'] = dict(labels=[16, 86, 142],eps_factor=0.2,L=50,wh=850,temp=293,obs='rate',pH=7.0,fasta=list(fasta_polyR30),ionic=0.02,path='./') 
    proteins.loc['polyU30'] = dict(labels=[16, 86, 142],eps_factor=0.2,L=50,wh=850,temp=300,obs='rate',pH=7.0,fasta=list(fasta_polyU30),ionic=0.02,path='./')
    return proteins

def genDCD(residues,name,prot,path,run_type,n_chains):
    """ 
    Generates coordinate and trajectory 
    in convenient formats for multiple chains
    """
    n_chains = int(n_chains)
    top = md.Topology()
    for _ in range(n_chains):
        chain = top.add_chain()
        for resname in prot.fasta:
            residue = top.add_residue(residues.loc[resname,'three'], chain)
            top.add_atom(residues.loc["p",'three'], element=md.element.phosphorus, residue=residue)
            top.add_atom(residues.loc[resname,'three'], element=md.element.nitrogen, residue=residue)
        for i in range(0,chain.n_atoms,2):
            if i==chain.n_atoms-2:
                top.add_bond(chain.atom(i),chain.atom(i+1))
            else:
                top.add_bond(chain.atom(i),chain.atom(i+1))
                top.add_bond(chain.atom(i),chain.atom(i+2))
    traj = md.load_dcd(path+"/{:s}.dcd".format(name), top)
    traj.center_coordinates()
    traj.xyz *= 1
    traj.unitcell_lengths *= 1
    traj.xyz += traj.unitcell_lengths[0,0]/2
    traj = traj.image_molecules(inplace=False, anchor_molecules=[set(traj.top.chain(0).atoms)],
           other_molecules=[set(traj.top.chain(1).atoms)], make_whole=False)
    traj[:].save_dcd(path+"/{:s}.dcd".format(run_type))
    traj[0].save_pdb(path+"/{:s}.pdb".format(run_type))

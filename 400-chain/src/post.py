
from analyse import *
import os
import subprocess
from jinja2 import Template
from rdf_calc import *
import argparse    

parser = argparse.ArgumentParser(description='fix trajctory and calculate rdf') 
parser.add_argument('arg1', help='name')  
parser.add_argument('arg2', help='n chains')  
args = parser.parse_args() 
name = args.arg1
n_chains = args.arg2

# protein template produced
proteins = initProteins()
proteins.to_pickle('proteins.pkl')
residues = pd.read_csv('../../residues.csv').set_index('one')
path = './'

# fix trajectories
genDCD(residues,name,proteins.loc[name],path,'fix',n_chains)

# calculate rdf
individual_rdfs(name,path,0.1,proteins,residues,n_chains)
concatenated_rdf(name,path,0.1,proteins,residues,n_chains)


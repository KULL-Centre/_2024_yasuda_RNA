import numpy as np
import os
from calvados.cfg import Config, Job, Components
import subprocess

cwd = os.getcwd()


config = Config(
  # GENERAL
  sysname = 'FUSRGG3', # name of simulation system
  box = [20,20,150], # nm
  temp = 293, # K
  ionic = 0.15, # molar
  pH = 7.5,
  topol = 'random',
  friction_coeff = 0.01,


  # RUNTIME SETTINGS
  wfreq = 1000000, # dcd writing frequency, 1 = 10fs
  steps = 1200000000, # number of simulation steps
  runtime = 0, # overwrites 'steps' keyword if > 0
  platform = 'CUDA', # 'CUDA'
  restart = 'checkpoint',
  verbose = True,

  # EQULIBRIATION SETTINGS
  slab_eq = True,
  k_eq = 0.01,
  steps_eq = 10000000,
)


# PATH
path = f'{cwd}/{config.config["sysname"]}'
subprocess.run(f'mkdir -p {path}',shell=True)

config.write(path,name='config.yaml')

components = Components(
  restraint = False, # apply restraints
  charge_termini = 'both', # charge N or C or both
  fresidues = f'{cwd}/residues.csv', # residue definitions
  ffasta = f'{cwd}/fastabib.fasta',

  # RNA settings
#  rna_kb1 = 1400.0,
#  rna_kb2 = 2200.0,
#  rna_ka = 4.20,
#  rna_pa = 3.14,
#  rna_nb_sigma = 0.4,
#  rna_nb_scale = 15,
#  rna_nb_cutoff = 2.0
)

components.add(name='FUSRGG3',  molecule_type='protein', nmol=500)
components.add(name='polyU40',  molecule_type='protein', nmol=100, charge_termini=False, kb=1000)

components.write(path,name='components.yaml')


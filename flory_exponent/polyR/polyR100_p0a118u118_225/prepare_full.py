
from os import getcwd
from calvados.cfg import preprocess

cwd = getcwd()

#########################################
## FULL SCRIPT ##
#########################################
fconfig = 'config.yaml'
ffasta = f'{cwd}/fastabib.fasta' # no fasta needed if pdb provided
batch_sys = 'PBS'
envname = 'calvados' # conda environment
fbash = '/home/people/sorbul/.bashrc'

# general simulation settings
box = [30., 30., 30.]
name = 'polyR100' #'polyA30'
temp = 293
ionic = .22500  #0.15
pH = 7.0

runtime_settings = {
    'wfreq' : 100000, # dcd writing frequency
    'steps' : 100000000, # number of simulation steps
    'runtime' : 0, # hours (overwrites steps if runtime is > 0)
    'pfname' : 'CPU', # 'CUDA',
    'threads' : 1, # number of threads for job script and openmm object
    'restart' : 'checkpoint', # 'pdb', None
    'frestart' : 'restart.chk', # restart file
}

general_settings = {
    'topol' : 'center',
    'cutoff_lj' : 2.0, # set the cutoff for the nonionic interactions
    'eps_factor' : 0.2, # eps_lj = eps_factor * 4.184
    'slab_eq' : False, # slab equilibration flag (restrain towards center)
    'k_eq' : 0.01, # kJ/mol*nm for linear restraints towards box center in z
    'res_csv' : '/home/yasuda/rna/exponent_factor/polyR30/polyR100_p0a118u118_225/residues.csv',
    'dt'      : 0.01 # time for one md step in ps 
}

# RNA parameters (single simulation)
rna_settings = {
    'rna' : {
        'polyR100': 1,
        },
    'rna_kb1' : 1400.0,          # bond force constant P-P bead, kJ/mol/nm^2, fitting 1000
    'rna_kb2' : 2200.0,          # bond force constant P-B bead
    'rna_ka'  : 4.2,             # angle force constant, kJ/mol 8
    'rna_pa'  : 3.14159265359,   # equilibrium angle phase, radian
    'rna_nb_sigma' : 0.4,        # LJ sigma for neighboring bases 
    'rna_nb_cutoff': 2.0         # LJ cutoff for neighboring bases
}


#########################################

loc = locals()
preprocess(loc)

#########################################

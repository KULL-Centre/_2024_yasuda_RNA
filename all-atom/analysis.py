import numpy as np
from scipy.optimize import curve_fit
import mdtraj as md

"""
UCAAUC
"""
t = md.load("UCAAUC/FFLJbb.dcd", top="UCAAUC/UCAAUC_AFORM.pdb")
top = t.topology

#P-P distance
p_ind = top.select("name P")
pairs = []
for i in range(1,len(p_ind)-2):
    pairs.append((p_ind[i], p_ind[i+1]))
pairs = np.array(pairs)
ppd = md.compute_distances(t, pairs)
np.save("UCAAUC_bond_pp.npy",ppd)

#P-base distance
pairs = []
pairs.append([top.select("name P and resid 1")[0], top.select("name N1 and resid 1")[0]])
pairs.append([top.select("name P and resid 2")[0], top.select("name N9 and resid 2")[0]])
pairs.append([top.select("name P and resid 3")[0], top.select("name N9 and resid 3")[0]])
pairs.append([top.select("name P and resid 4")[0], top.select("name N1 and resid 4")[0]])
#pairs.append([top.select("name P and resid 5")[0], top.select("name N1 and resid 5")[0]])
pairs = np.array(pairs)
ppd = md.compute_distances(t, pairs)
np.save("UCAAUC_bond_pbase.npy",ppd)

#p-p-p angle
pairs = []
pairs.append([top.select("name P and resid 1")[0],
              top.select("name P and resid 2")[0],
              top.select("name P and resid 3")[0]])
pairs.append([top.select("name P and resid 2")[0],
              top.select("name P and resid 3")[0],
              top.select("name P and resid 4")[0]])
#pairs.append([top.select("name P and resid 3")[0],
#              top.select("name P and resid 4")[0],
#              top.select("name P and resid 5")[0]])
pairs = np.array(pairs)
ppd = md.compute_angles(t, pairs)
np.save("UCAAUC_angle_ppp.npy",ppd)


#base-P-base-P dihedral
pairs = []
pairs.append([top.select("name N1 and resid 1")[0],
                         top.select("name P and resid 1")[0],
                         top.select("name P and resid 2")[0],
                         top.select("name N9 and resid 2")[0]])
pairs.append([top.select("name N9 and resid 2")[0],
                         top.select("name P and resid 2")[0],
                         top.select("name P and resid 3")[0],
                         top.select("name N9 and resid 3")[0]])
pairs.append([top.select("name N9 and resid 3")[0],
                         top.select("name P and resid 3")[0],
                         top.select("name P and resid 4")[0],
                         top.select("name N1 and resid 4")[0]])             
#pairs.append([top.select("name N1 and resid 4")[0],
#                         top.select("name P and resid 4")[0],
#                         top.select("name P and resid 5")[0],
#                         top.select("name N1 and resid 5")[0]])
pairs = np.array(pairs)
ppd = md.compute_dihedrals(t, pairs)
np.save("UCAAUC_dihedral_baseppbase.npy",ppd)

# base-base distance
pairs = []
pairs.append([top.select("name N1 and resid 1")[0],
              top.select("name N9 and resid 2")[0]])
pairs.append([top.select("name N9 and resid 2")[0],
              top.select("name N9 and resid 3")[0]])
pairs.append([top.select("name N9 and resid 3")[0],
              top.select("name N1 and resid 4")[0]])             
#pairs.append([top.select("name N1 and resid 4")[0],
#                         top.select("name N1 and resid 5")[0]])
ppd = md.compute_distances(t, pairs)
np.save("UCAAUC_bond_basebase.npy",ppd)

"""
rA30 model
"""

t = md.load("./HCG_models/rA30_20models.pdb") # poly rA30
top = t.topology

# p-p distance
p_ind = top.select("name P")
pairs = []
for i in range(1,len(p_ind)-2):
    pairs.append((p_ind[i], p_ind[i+1]))
pairs = np.array(pairs)
ppd = md.compute_distances(t, pairs)
np.save("rA30_bond_pp.npy",ppd)

#p-base distance
pairs = []
for i in range(1,29):
    pairs.append([top.select("name P and resid {:d}".format(i))[0], 
                  top.select("name N9 and resid {:d}".format(i))[0]])
pairs = np.array(pairs)
ppd = md.compute_distances(t, pairs)
np.save("rA30_bond_pbase.npy",ppd)

pairs = []
for i in range(1,28):
    pairs.append([top.select("name N9 and resid {:d}".format(i))[0],
                  top.select("name P and resid {:d}".format(i))[0],
                  top.select("name P and resid {:d}".format(i+1))[0],
                  top.select("name N9 and resid {:d}".format(i+1))[0]])
pairs = np.array(pairs)
ppd = md.compute_dihedrals(t, pairs)
np.save("rA30_dihedral_baseppbase.npy",ppd)

pairs = []
for i in range(1,27):
    pairs.append([top.select("name P and resid {:d}".format(i))[0],
                  top.select("name P and resid {:d}".format(i+1))[0],
                  top.select("name P and resid {:d}".format(i+2))[0]])
pairs = np.array(pairs)
ppd = md.compute_angles(t, pairs)
np.save("rA30_angle_ppp.npy",ppd)

#base-base distance
pairs = []
for i in range(1,28):
    pairs.append([top.select("name N9 and resid {:d}".format(i))[0],
                  top.select("name N9 and resid {:d}".format(i+1))[0]])
pairs = np.array(pairs)
ppd = md.compute_distances(t, pairs)
np.save("rA30_bond_basebase.npy",ppd)


"""
rU30 model
"""

t = md.load("./HCG_models/rU30_20models.pdb") # poly rA30
top = t.topology

# p-p distance
p_ind = top.select("name P")
pairs = []
for i in range(1,len(p_ind)-2):
    pairs.append((p_ind[i], p_ind[i+1]))
pairs = np.array(pairs)
ppd = md.compute_distances(t, pairs)
np.save("rU30_bond_pp.npy",ppd)

#p-base distance
pairs = []
for i in range(1,29):
    pairs.append([top.select("name P and resid {:d}".format(i))[0],
                  top.select("name N1 and resid {:d}".format(i))[0]])
pairs = np.array(pairs)
ppd = md.compute_distances(t, pairs)
np.save("rU30_bond_pbase.npy",ppd)

pairs = []
for i in range(1,28):
    pairs.append([top.select("name N1 and resid {:d}".format(i))[0],
                  top.select("name P and resid {:d}".format(i))[0],
                  top.select("name P and resid {:d}".format(i+1))[0],
                  top.select("name N1 and resid {:d}".format(i+1))[0]])
pairs = np.array(pairs)
ppd = md.compute_dihedrals(t, pairs)
np.save("rU30_dihedral_baseppbase.npy",ppd)

pairs = []
for i in range(1,27):
    pairs.append([top.select("name P and resid {:d}".format(i))[0],
                  top.select("name P and resid {:d}".format(i+1))[0],
                  top.select("name P and resid {:d}".format(i+2))[0]])
pairs = np.array(pairs)
ppd = md.compute_angles(t, pairs)
np.save("rU30_angle_ppp.npy",ppd)

#base-base distance
pairs = []
for i in range(1,28):
    pairs.append([top.select("name N1 and resid {:d}".format(i))[0],
                  top.select("name N1 and resid {:d}".format(i+1))[0]])
pairs = np.array(pairs)
ppd = md.compute_distances(t, pairs)
np.save("rU30_bond_basebase.npy",ppd)





























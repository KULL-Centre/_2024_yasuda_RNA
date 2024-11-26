import mdtraj as md
import numpy as np
import pandas as pd
def read_exp(_name="polyU30",fn='./exp/'):
    f = fn+"pollack2020.csv"
    df = pd.read_csv(f)
    exp = df[df['name']==_name]
    return exp

def calcRg(t,masses):
    # calculate the center of mass
    cm = np.sum(t.xyz*masses[np.newaxis,:,np.newaxis],axis=1)/masses.sum()
    # calculate residue-cm distances
    si = np.linalg.norm(t.xyz - cm[:,np.newaxis,:],axis=2)
    # calculate rg
    rgarray = np.sqrt(np.sum(si**2*masses,axis=1)/masses.sum())
    return rgarray


def read_sim(_name,_p,_a,_u,skip=500): #100ps*500=50 ns
    T=293
    df = []
    for i in ions:
        key=_name+"_p"+str(_p)+"a"+str(_a)+"u"+str(_u)+"_"+str(i)
        trj = md.load(dir_+"/{0}/{1}/{2}/{1}.dcd".format(key,_name,T), 
                      top=dir_+"/{0}/{1}/{2}/top.pdb".format(key,_name,T))  
        rg = calcRg(trj, masses=get_mass(_name))[skip:]
        mean,std = np.mean(rg),np.std(rg)
        df.append([i,mean,std])
    df = pd.DataFrame(df, columns=["conc","mean", "std"])
    return df

def get_mass(_name,n=30,mp=194.1,ma=134.1,mu=111.1,mr=126.3):
    mass = []
    if _name=="polyA30":
        for i in range(n):
            mass.append(mp)
            mass.append(ma)
    elif _name=="polyU30": 
        for i in range(n):
            mass.append(mp)
            mass.append(mu)
    elif _name[:5]=="polyU":
        for i in range(n):
            mass.append(mp)
            mass.append(mu)
    elif _name[:5]=="polyR":
        for i in range(n):
            mass.append(mp)
            mass.append(mr)
    return np.array(mass)

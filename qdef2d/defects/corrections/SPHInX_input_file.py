import os
import json
import numpy as np
import argparse
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pymatgen
from pymatgen.core.lattice import Lattice


def Ang_to_bohr(x):
    
    ## convert from Angstroms to bohr
    return (x * 1.88973)


def structure_grp(lattice):      
    
    s = 'structure { \n'
    s += '    cell = '+ str(Ang_to_bohr(lattice.matrix).tolist()) + '; \n'
    s += '} \n\n'

    return s


def slab_grp(slabmin,slabmax,eps):

    s = 'slab { \n'
    s += '    fromZ = {}; \n'.format(Ang_to_bohr(slabmin))
    s += '    toZ = {}; \n'.format(Ang_to_bohr(slabmax))
    s += '    epsilon = {}; \n'.format(np.trace(eps)/3)
    s += '} \n\n'  

    return s


def charge_grp(posZ,q):

    s = 'charge { \n'
    s += '    posZ = {}; \n'.format(Ang_to_bohr(posZ))
    s += '    Q = {:+}; \n'.format(-q)
    s += '} \n\n'  

    return s


def isolated_grp(slabmin,slabmax):
    
    if slabmin <= 3.0:
        raise ValueError ("insufficient vacuum to define isolated group")

    s = 'isolated { \n'
    s += '    fromZ = {}; \n'.format(Ang_to_bohr(slabmin-3))
    s += '    toZ = {}; \n'.format(Ang_to_bohr(slabmax+3))
    s += '} \n\n'  

    return s

    
def generate(eps,slab_d):
    
    """
    Generate SPHInX input file.
    
    eps (float): averaged dielectric constant
    slab_d (float): corresponding slab thickness (Angstroms)
    
    """

    ## the main script should have put us in the appropriate subdirectory
    dir_sub = os.getcwd()
              
    with open(os.path.join(dir_sub,"defectproperty.json"), 'r') as file:
        defprop = json.loads(file.read())
    lattice = Lattice.from_dict(defprop["lattice"])
    
    ## STRUCTURE GROUP
    s = structure_grp(lattice)
    
    ## SLAB GROUP
    ## assume slab is centered in the cell vertically
    slabmin = (lattice.c-slab_d)/2
    slabmax = (lattice.c+slab_d)/2
    s += slab_grp(slabmin,slabmax,eps=eps*np.eye(3))
    
    ## CHARGE GROUP
    posZ = np.mean([def_site[2] for def_site in defprop["defect_site"]])
    s += charge_grp(posZ*lattice.c,defprop["charge"])
    
    ## ISOLATED GROUP
    s += isolated_grp(slabmin,slabmax)

    if not os.path.exists(os.path.join(dir_sub,"correction")):
            os.makedirs(os.path.join(dir_sub,"correction"))
    with open(os.path.join(dir_sub,"correction","system.sx"),'w') as f:
        f.write(s)
        
        
if __name__ == '__main__':
    

    ## this script can also be run directly from the command line
    parser = argparse.ArgumentParser(description='Generate SPHInX input file')
    parser.add_argument('eps',type=float,help='averaged dielectric constant')
    parser.add_argument('slab_d',type=float,help='corresponding slab thickness (Angstroms)')
     
    ## read in the above arguments from command line
    args = parser.parse_args()
    
    generate(args.eps, args.slab_d)
    
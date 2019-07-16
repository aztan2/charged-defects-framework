import os
import json
import numpy as np
import argparse
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pymatgen
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.io.vasp.inputs import Poscar


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


if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description='Generate INCAR')
    parser.add_argument('eps',type=float,help='averaged dielectric constant')
    parser.add_argument('slab_d',type=float,help='corresponding slab thickness (Angstroms)')

      
    ## read in the above arguments from command line
    args = parser.parse_args()
    

    ## the bash script already put us in the appropriate subdirectory
    dir_sub = os.getcwd()

#    lattice = Poscar.from_file(folder1+"POSCAR").structure.lattice                
    with open(os.path.join(dir_sub,"defectproperty.json"), 'r') as file:
        defprop = json.loads(file.read())
    lattice = Lattice.from_dict(defprop["lattice"])
    
    ## STRUCTURE GROUP
    s = structure_grp(lattice)
    
    ## SLAB GROUP
    ## assume slab is centered in the cell vertically
    slabmin = (lattice.c-args.slab_d)/2
    slabmax = (lattice.c+args.slab_d)/2
    s += slab_grp(slabmin,slabmax,eps=args.eps*np.eye(3))
    
    ## CHARGE GROUP
    posZ = np.mean([def_site[2] for def_site in defprop["defect_site"]])
    s += charge_grp(posZ*lattice.c,defprop["charge"])
    
    ## ISOLATED GROUP
    s += isolated_grp(slabmin,slabmax)

    if not os.path.exists(os.path.join(dir_sub,"correction")):
            os.makedirs(os.path.join(dir_sub,"correction"))
    with open(os.path.join(dir_sub,"correction","system.sx"),'w') as f:
        f.write(s)
                        
    
    
##    folder = "../mp-2815_MoS2/test/GGA/mag/"
#    folder = "Y:/MoS2/monolayer_Nbsub/SCAN_vdW/mag/"
#
#    qs = [-2]
#    cells = [(4,4),(5,5)]
#    vacs = [15,20]
#
#
##    eps_ave = 17.18  ## averaged "isotropic" dielectric constant (PBE)
#    eps_ave = 16.27  ## averaged "isotropic" dielectric constant (SCAN w/o ionic)
#    slab_d = 5.40  ## corresponding slab thickness in Angstroms 
##    atomrad = 1.14  ## Sulphur atomic radius ~ 100pm = 1 Angstrom
#    
#    
#    for q in qs:
#        for cell in cells:
#            for vac in vacs:
#                
#                folder1 = folder+"charge_%d/%dx%dx1/vac_%d/"%(q,cell[0],cell[1],vac)
##                folder1 = folder+"charge_%d/%dx%dx1/vac_%d/ismear1/"%(q,cell[0],cell[1],vac)
#
##                lattice = Poscar.from_file(folder1+"POSCAR").structure.lattice                
#                with open(folder1+"defectproperty.json", 'r') as file:
#                    defprop = json.loads(file.read())
#                lattice = Lattice.from_dict(defprop["lattice"])
#                
#                ## STRUCTURE GROUP
#                s = structure_grp(lattice)
#                
#                ## SLAB GROUP
#                ## assume slab is centered in the cell vertically
#                slabmin = (lattice.c-slab_d)/2
#                slabmax = (lattice.c+slab_d)/2
#                s += slab_grp(slabmin,slabmax,eps=eps_ave*np.eye(3))
#                
#                ## CHARGE GROUP
#                posZ = np.mean([def_site[2] for def_site in defprop["defect_site"]])
#                s += charge_grp(posZ*lattice.c,q)
#                
#                ## ISOLATED GROUP
#                s += isolated_grp(slabmin,slabmax)
#    
##                folder2 = folder + "../../../corr_sxdefectalign2d/GGA/varyeps/d3.12/%dx%dx1/vac_%d/charge_%d/"%(cell[0],cell[1],vac,q)
#                folder2 = "Y:/MoS2/corrections/sxdefectalign2d/monolayer_Nbsub/SCAN_vdW/mag/%dx%dx1/vac_%d/charge_%d/"%(cell[0],cell[1],vac,q)
##                folder2 = folder1
#                if not os.path.exists(folder2):
#                    os.makedirs(folder2)
#                with open(folder2+"system.sx",'w') as f:
#                    f.write(s)
#                    
                    
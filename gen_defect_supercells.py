import os
import json
import argparse
import numpy as np
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pymatgen
from pymatgen.io.vasp.inputs import Poscar


class Defect():
    
    def __init__(self, structure, supercell, vacuum, charge):
        
        self.structure = structure
        self.supercell = supercell
        self.vacuum = vacuum
        self.charge = charge
        self.defect_type = []
        self.defect_site = []
            
        
    def add_defect_info(self, defect_type, defect_site):
        
        self.defect_type.append(defect_type)
        self.defect_site.append(defect_site)
            
        return
    

    def remove_atom(self):

        for def_type,def_site in zip(self.defect_type,self.defect_site): 
            if def_type[0] == "v":
                if def_site in self.structure:
                    defect_index = self.structure.index(def_site)
                    self.structure.remove_sites([defect_index])
                else:
                    raise ValueError ("site does not exist to create a vacancy")
        
        return


    def replace_atom(self):
     
        for def_type,def_site in zip(self.defect_type,self.defect_site):
            print (def_type,def_site)
            if def_type[0] == "s":
                if def_site in self.structure:
                    defect_index = self.structure.index(def_site)
                    self.structure[defect_index] = def_type.split("_")[-1]
                else:
                    raise ValueError ("site does not exist to create a substitutional")
        
        return
    
    
    def as_dict(self):
        
        d = {"charge": self.charge,
             "defect_type": self.defect_type,
             "defect_site": [site.as_dict()["abc"] for site in self.defect_site],
             "lattice": self.structure.lattice.as_dict(),
             "supercell": self.supercell,
             "vacuum": self.vacuum
             }
                    
        return d

    
if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description='Generate defect supercells.')
    parser.add_argument('dir_poscars',help='path to directory with unitcell POSCARs')
    parser.add_argument('q',type=int,help='charge')
    parser.add_argument('supercell',help='supercell size')
    parser.add_argument('vacuum',type=int,help='vacuum spacing')
    parser.add_argument('--write_bulkref',help='to write bulkref?',default=False,action='store_true')
      
    ## read in the above arguments from command line
    args = parser.parse_args()
    
    
    ## the bash script already put us in the appropriate subdirectory
    dir_sub = os.getcwd()
    
    ## read in corresponding unit cell POSCAR
    poscar = Poscar.from_file(os.path.join(args.dir_poscars,"POSCAR_vac_%d"%args.vacuum),
                              check_for_POTCAR=False, read_velocities=False)       
    
    ## make undefected supercell
    structure = poscar.structure.copy()
    nvecs = [int(n) for n in args.supercell.split('x')]
    structure.make_supercell([nvecs[0],nvecs[1],nvecs[2]])
    structure_bulk = structure.copy()
    
    ## initialize defect object
    defect = Defect(structure.copy(),supercell=nvecs,vacuum=args.vacuum,charge=args.q)
        
    ## create vacancy defect
    ## THIS PART IS STILL HARDCODED !!
#    defect_site = structure_bulk[structure_bulk.num_sites-1]
#    defect.add_defect_info(defect_type="vac_S",defect_site=defect_site)
    defect_site = structure_bulk[int(structure_bulk.num_sites/3-1)]
    defect.add_defect_info(defect_type="vac_W",defect_site=defect_site)
    ## create second vacancy defect (divacancy)
#    defect_site = structure_bulk[int(2*structure_bulk.num_sites/3-1)]
#    defect_site = structure_bulk[structure_bulk.num_sites-2]
#    defect.add_defect_info(defect_type="vac_Se",defect_site=defect_site)
    defect.remove_atom()
        
    ## create substitutional defect
#    defect_site = structure_bulk[int(structure_bulk.num_sites/3-1)]
#    defect.add_defect_info(defect_type="sub_Nb",defect_site=defect_site)
#    defect_site = structure_bulk[int(structure_bulk.num_sites-1)]
#    defect.add_defect_info(defect_type="sub_O",defect_site=defect_site)
#    defect.replace_atom()

    ## write POSCAR
    Poscar.write_file(Poscar(defect.structure.get_sorted_structure()),os.path.join(dir_sub,"POSCAR"))
    
    ## write bulkref POSCAR
    if args.write_bulkref:
        if args.q == 0:
            if not os.path.exists(os.path.join(dir_sub,"bulkref")):
                os.makedirs(os.path.join(dir_sub,"bulkref"))
            Poscar.write_file(Poscar(structure_bulk),os.path.join(dir_sub,"bulkref","POSCAR"))
      
    ## write json file
    with open(os.path.join(dir_sub,"defectproperty.json"), 'w') as file:
         file.write(json.dumps(defect.as_dict(),indent=4)) # use `json.loads` to do the reverse
    
    
##    dir_main = "mp-2815_MoS2/test/GGA/mag/"
#    dir_main = "Y:/MoS2/monolayer_Nbsub/GGA/mag/"
#    if not os.path.exists(dir_main):
#        os.makedirs(dir_main)
#
#
#    for vac in [15,20]:
#        ## read in corresponding unit cell POSCAR
#        poscar = Poscar.from_file("mp-2815_MoS2/monolayer/GGA/mag/vac_%d/POSCAR"%vac)
##        poscar = Poscar.from_file("mp-2815_MoS2/monolayer/GGA/mag/vac_%d/POSCAR_orth"%vac)
#            
#        for nvecs in [[4,4,1],[5,5,1]]:   
#            ## make undefected supercell
#            structure = poscar.structure.copy()
#            structure.make_supercell([nvecs[0],nvecs[1],nvecs[2]])
#            structure_bulk = structure.copy()
#                
#            for q in [-2]:
#                ## create the subdirectory
#                dir_sub = dir_main+"charge_%d/%dx%dx%d/vac_%d/"%(q,nvecs[0],nvecs[1],nvecs[2],vac)
#                if not os.path.exists(dir_sub):
#                    os.makedirs(dir_sub)
#
#                ## initialize defect object
#                defect = Defect(structure.copy(),supercell=nvecs,vacuum=vac,charge=q)
#                
#                ## create S vacancy defect
##                defect_site = structure_bulk[structure_bulk.num_sites-1]
##                defect.add_defect_info(defect_type="vac_S",defect_site=defect_site)
##                defect.remove_atom()
#                
#                ## create substitutional defect
#                defect_site = structure_bulk[int(structure_bulk.num_sites/3-1)]
#                defect.add_defect_info(defect_type="sub_Nb",defect_site=defect_site)
#                defect.replace_atom(new_species="Nb")
#    
#                ## write POSCAR
#                Poscar.write_file(Poscar(defect.structure.get_sorted_structure()),dir_sub+"POSCAR")
#                if q == 0:
#                    if not os.path.exists(dir_sub+"/bulkref"):
#                        os.makedirs(dir_sub+"/bulkref")
#                    Poscar.write_file(Poscar(structure_bulk),dir_sub+"bulkref/POSCAR")
#  
#                ## write json file
#                with open(dir_sub+"defectproperty.json", 'w') as file:
#                     file.write(json.dumps(defect.as_dict(),indent=4)) # use `json.loads` to do the reverse


    ## BULK  
#    dir_main = "mp-2815_MoS2/test/GGA/mag/charge_0/"
#    if not os.path.exists(dir_main):
#        os.makedirs(dir_main)
#    
#    for nvecs in [[3,3,1]]:
#        dir_sub = dir_main+"%dx%dx%d/"%(nvecs[0],nvecs[1],nvecs[2])
#            
#        ## read in corresponding unit cell POSCAR
#        poscar = Poscar.from_file("mp-2815_MoS2/bulk/GGA/expt_c/POSCAR")
#        
#        ## make undefected supercell
#        structure = poscar.structure
#        structure.make_supercell([nvecs[0],nvecs[1],nvecs[2]])
#        structure_bulk = structure.copy()
#        if not os.path.exists(dir_sub+"/bulkref"):
#            os.makedirs(dir_sub+"/bulkref")
#        Poscar.write_file(Poscar(structure_bulk),dir_sub+"bulkref/POSCAR")
#        
#        ## create S vacancy defect
#        defect_site = structure_bulk[structure.num_sites-1]
#        defect = Defect(structure,defect_type="vac_S",defect_site=defect_site,charge=0)
#        defect.remove_atom()
#        
#        ## write defect POSCARS            
#        if not os.path.exists(dir_sub):
#            os.makedirs(dir_sub)
#        Poscar.write_file(Poscar(defect.structure.get_sorted_structure()),dir_sub+"POSCAR")
#
#        ## write json file
#        with open(dir_sub+"defectproperty.json", 'w') as file:
#             file.write(json.dumps(defect.as_dict(),indent=4)) # use `json.loads` to do the reverse

    

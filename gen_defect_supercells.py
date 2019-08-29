import os
import json
import argparse
import itertools
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
            

    def get_defect_site(self, initdef):

        ## figuring out min and max+1 indices for each species        
        syms = [site.specie.symbol for site in self.structure]
        elems = [el[0] for el in itertools.groupby(syms)]
        elem_max = np.cumsum([len(tuple(el[1])) for el in itertools.groupby(syms)])
        elem_min = [0] + [i for i in elem_max[:-1]]
        elem_minmax = {el: [imin,imax] for el,imin,imax in zip(elems,elem_min,elem_max)}
        
        if initdef["type"][0] == "v" or initdef["type"][0] == "s":
            ## initdef["index"] is a relative site index
            ## if possible, figure out a less clunky way to index...
            if initdef["index"] < 0:
                defect_site = structure_bulk[int(elem_minmax[initdef["species"]][1]+initdef["index"])]
            else:
                defect_site = structure_bulk[int(elem_minmax[initdef["species"]][0]+initdef["index"])]
            
        return defect_site
    
        
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
    parser.add_argument('json_initdef',help='json file with details to initialize defect')
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
    
    
    ## read in defect details from initdefect json file
    with open(args.json_initdef, 'r') as file:
        initdef = json.loads(file.read())

    ## initialize defect object
    defect = Defect(structure.copy(),supercell=nvecs,vacuum=args.vacuum,charge=args.q)

    ## set the defect info (type, site, species) for each defect listed in the json file
    for d in initdef:
        defect_site = defect.get_defect_site(initdef[d])
        if initdef[d]["type"][0] == "v": 
            ## vacancy type defect
            defect.add_defect_info(defect_type=initdef[d]["type"]+"_"+initdef[d]["species"],
                                   defect_site=defect_site)
        if initdef[d]["type"][0] == "s": 
            ## substitutional type defect
            defect.add_defect_info(defect_type=initdef[d]["type"]+"_"+initdef[d]["species_new"],
                                   defect_site=defect_site)            
#        defect_site = structure_bulk[int(2*structure_bulk.num_sites/3-1)]

    ## create vacancy defect(s)           
    defect.remove_atom()
    ## create substitutional defect(s)
    defect.replace_atom()
        

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

    

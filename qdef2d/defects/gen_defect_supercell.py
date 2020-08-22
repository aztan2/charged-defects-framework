import sys
import os
import json
import argparse
import numpy as np
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pymatgen
from pymatgen.io.vasp.inputs import Poscar
from qdef2d.defects.core import Defect


def main(args):
    
    ## define a main function callable from another python script
    
    parser = argparse.ArgumentParser(description='Generate defect supercells.')
    parser.add_argument('dir_poscars',help='path to directory with unitcell POSCARs')
    parser.add_argument('json_initdef',help='json file with details to initialize defect')
    parser.add_argument('q',type=int,help='charge')
    parser.add_argument('supercell',help='supercell size')
    parser.add_argument('vacuum',type=int,help='vacuum spacing')
    parser.add_argument('--write_bulkref',help='to write bulkref?',
                        default=False,action='store_true')
      
    ## parse the given arguments
    args = parser.parse_args(args)

    ## this script should be called from within
    ## the appropriate charge/cell/vacuum subdirectory
    ## get the current working directory
    dir_sub = os.getcwd()
    
    
    ## read in corresponding unit cell POSCAR
    poscar = Poscar.from_file(os.path.join(args.dir_poscars,"POSCAR_vac_%d"%args.vacuum),
                              check_for_POTCAR=False, read_velocities=False)     
    
    ## make undefected supercell
    structure = poscar.structure.copy()
    nvecs = [int(n) for n in args.supercell.split('x')]
    structure.make_supercell([nvecs[0],nvecs[1],nvecs[2]])
    structure_bulk = structure.copy()
    
    if args.write_bulkref:
        ## write bulkref POSCAR
        Poscar.write_file(Poscar(structure_bulk),
                          os.path.join(dir_sub,"POSCAR"))
        
    else:
        ## initialize defect object
        defect = Defect(structure_bulk,structure.copy(),nvecs,args.vacuum,args.q)
        
        ## read in defect details from initdefect json file
        with open(args.json_initdef, 'r') as file:
            initdef = json.loads(file.read())
    
        ## set the defect info (type, site, species) for each defect listed in the json file
        for d in initdef:
            if defect.parse_initdef(initdef[d]):
                ## if index offsets are not specified in the json file, set them to be 0
                initdef[d]["index_offset_a"] = initdef[d].get("index_offset_a",0)
                initdef[d]["index_offset_b"] = initdef[d].get("index_offset_b",0)
                initdef[d]["index_offset_basis"] = initdef[d].get("index_offset_basis",0)
                print (initdef[d])
                defect_site = defect.get_defect_site(initdef[d])
                defect.add_defect_info(initdef[d],defect_site)  
    
        ## create vacancy defect(s)           
        defect.remove_atom()
        ## create substitutional defect(s)
        defect.replace_atom()
        ## create interstitial defect(s)
        defect.add_atom()     
    
        ## write defect POSCAR
        Poscar.write_file(Poscar(defect.structure.get_sorted_structure()),
                          os.path.join(dir_sub,"POSCAR"))
          
        ## write json file
        with open(os.path.join(dir_sub,"defectproperty.json"), 'w') as file:
             file.write(json.dumps(defect.as_dict(),indent=4))  
    
    
if __name__ == '__main__':

    main(sys.argv[1:])    

    

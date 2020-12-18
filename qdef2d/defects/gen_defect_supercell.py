import os
import errno
import json
import argparse
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pymatgen
from pymatgen.io.vasp.inputs import Poscar
from qdef2d.defects.core import Defect


def generate(dir_def_main, initdef_file, q, supercell, vacuum, bulkref=False):

    """ 
    Generate defect supercell.
    
    Parameters
    ----------
    dir_def_main (str): path to main defect directory containing unitcell POSCARs
    initdef_file (str): json file with details to initialize defect
    q (int): charge
    supercell (tuple of ints): supercell size as [n1,n2,n3]
    vacuum (int): vacuum spacing
    [optional] bulkref (bool): generate only bulk reference supercell? Default=False.
    
    """
    

    ## the directory from which this function was called
    ## should already be the appropriate charge/cell/vacuum subdirectory
    subdir_def = os.getcwd()
    
    
    ## read in corresponding unit cell POSCAR
    pos_file = os.path.join(dir_def_main,"POSCAR_vac_%d"%vacuum)
    if not os.path.exists(pos_file):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), pos_file)
    poscar = Poscar.from_file(pos_file,check_for_POTCAR=False,read_velocities=False)
    
    
    ## make undefected supercell
    structure = poscar.structure.copy()
#    nvecs = [int(n) for n in supercell.split('x')]
    structure.make_supercell(supercell)
    structure_bulk = structure.copy()
    
    
    if bulkref:
        ## write bulkref POSCAR
        Poscar.write_file(Poscar(structure_bulk),
                          os.path.join(subdir_def,"POSCAR"))
        
    else:
        ## initialize defect object
        defect = Defect(structure_bulk,structure.copy(),supercell,vacuum,q)
        
        ## read in defect details from initdefect json file
        id_file = os.path.join(dir_def_main,initdef_file)
        if not os.path.exists(id_file):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), id_file)
        with open(id_file, 'r') as file:
            initdef = json.loads(file.read())
    
        ## for each defect listed in the json file
        for d in initdef:
            ## parse details from initdef, 
            ## checking for required fields and correct entry types
            initdef_d = defect.parse_initdef(initdef[d])
#            print (initdef_d)
            ## determine the defect site (PeriodicSite object)
            defect_site = defect.get_defect_site(initdef_d)
            ## set the defect info (type, species, site) 
            defect.add_defect_info(initdef_d,defect_site)  
    
        ## create vacancy defect(s)           
        defect.remove_atom()
        ## create substitutional defect(s)
        defect.replace_atom()
        ## create interstitial defect(s)
        defect.add_atom()     
    
        ## write defect POSCAR
        Poscar.write_file(Poscar(defect.structure.get_sorted_structure()),
                          os.path.join(subdir_def,"POSCAR"))
          
        ## write defectproperty.json file (summary of defect created in this supercell)
        with open(os.path.join(subdir_def,"defectproperty.json"), 'w') as file:
             file.write(json.dumps(defect.as_dict(),indent=4))  
    
    
if __name__ == '__main__': 
    
    
    ## this script can also be run directly from the command line
    parser = argparse.ArgumentParser(description='Generate defect supercell.')
    parser.add_argument('dir_def_main',help='path to main defect directory containing unitcell POSCARs')
    parser.add_argument('initdef_file',help='json file with details to initialize defect')
    parser.add_argument('q',type=int,help='charge')
    parser.add_argument('supercell',help='supercell size, as n1xn2xn3')
    parser.add_argument('vacuum',type=int,help='vacuum spacing')
    parser.add_argument('--bulkref',help='generate only bulk reference supercell?',
                        default=False,action='store_true')
      
    ## parse the given arguments
    args = parser.parse_args()


    generate(args.dir_def_main, args.initdef_file, 
             args.q, [int(n) for n in args.supercell.split('x')], args.vacuum, 
             args.bulkref)  

import sys
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
from pymatgen.core.sites import Site,PeriodicSite


class Defect():
    
    def __init__(self, structure_bulk, structure, supercell, vacuum, charge):
        
        self.structure_bulk = structure_bulk
        self.structure = structure
        self.supercell = supercell
        self.vacuum = vacuum
        self.charge = charge
        self.defect_type = []
        self.defect_site = []

   
    def get_site_ind(self, ind, sp, offset_n1, offset_n2, offset_n1n2):

        ## figuring out min and max+1 indices for each species        
        syms = [site.specie.symbol for site in self.structure]
        elems = [el[0] for el in itertools.groupby(syms)]
        elem_max = np.cumsum([len(tuple(el[1])) for el in itertools.groupby(syms)])
        elem_min = [0] + [i for i in elem_max[:-1]]
        elem_minmax = {el: [imin,imax] for el,imin,imax in zip(elems,elem_min,elem_max)}

        ## initdef["index"] is a relative site index
        ## if possible, figure out a less clunky way to index...
        siteind_new = (ind + offset_n1n2*self.supercell[0]*self.supercell[1]
                           + offset_n1*self.supercell[1] + offset_n2)
        if ind < 0:
            siteind_new = int(elem_minmax[sp][1]+siteind_new)
        else:
            siteind_new = int(elem_minmax[sp][0]+siteind_new)   
        
        return siteind_new
    

    def get_defect_site(self, initdef):

        siteinds = []
        
        ## this is relatively simple for vacancies and substitutionals
        ## which by definition are a pre-existing site in the structure        
        if initdef["type"][0] == "v" or initdef["type"][0] == "s":
            siteind_new = self.get_site_ind(initdef["index"],
                                            initdef["species"],
                                            initdef["index_offset_n1"],
                                            initdef["index_offset_n2"],
                                            initdef["index_offset_n1n2"])
            siteinds.append(siteind_new)
            defect_site = self.structure_bulk[siteind_new]
        
        ## for interstitials/adatoms, the defect site coords are defined
        ## by the center (average) of the user-provided list of sites
        if initdef["type"][0] == "i" or initdef["type"][0] == "a":
            if len(initdef["index"]) != len(initdef["species"]):
                raise ValueError ("inconsistency between index and species lists")
            for offset in ["index_offset_n1","index_offset_n2","index_offset_n1n2"]:
                if initdef[offset] == 0:
                    initdef[offset] = [0]*len(initdef["index"])

            ## get absolute site indices
            siteind_new = [self.get_site_ind(ind,sp,offset_n1,offset_n2,offset_n1n2)
                            for ind,sp,offset_n1,offset_n2,offset_n1n2
                                in zip(initdef["index"],initdef["species"],
                                       initdef["index_offset_n1"],
                                       initdef["index_offset_n2"],
                                       initdef["index_offset_n1n2"])]
            siteinds.append(siteind_new)
            
            ## get the averaged position
            coords_ref = self.structure_bulk[siteind_new[0]].frac_coords
            defect_coords = coords_ref
            for i in siteind_new[1:]:
                coords = self.structure_bulk[i].frac_coords
                for j in [0,1,2]:
                    if abs(coords[j] - coords_ref[j]) > 0.5:  ## dealing with pbc 
                        if coords[j] > coords_ref[j]: coords[j] -= 1.0
                        else: coords[j] += 1.0
                defect_coords += coords
            defect_coords = defect_coords/len(siteind_new)

            ## we may want to shift the z position, e.g. in the case of an adatom            
            if initdef.get("shift_z") != None:
                defect_coords[2] += float(initdef["shift_z"])/self.structure.lattice.c
            
            ## create the defect_site as a PeriodicSite object
            defect_site = PeriodicSite(initdef["species_new"],defect_coords,
                                       lattice=self.structure.lattice,
                                       coords_are_cartesian=False)
        
        return defect_site, siteinds


    def add_defect_info(self, initdef, defect_site):
        
        if initdef["type"][0] == "v":
            ## vacancy type defect
            self.defect_type.append(initdef["type"]+"_"+initdef["species"])
            self.defect_site.append(defect_site)
        if initdef["type"][0] in ["s", "i", "a"]: 
            ## substitutional/interstitial/adatom type defect
            self.defect_type.append(initdef["type"]+"_"+initdef["species_new"])
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
    
    
    def add_atom(self):
        
        for def_type,def_site in zip(self.defect_type,self.defect_site):
            if def_type[0] == "i" or def_type[0] == "a":
                if def_site.frac_coords.tolist() in self.structure.frac_coords.tolist():
                    raise ValueError ("site already exists; can't add an atom here")
                else:
                    self.structure.append(def_site.specie,def_site.frac_coords)
                
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


def main(args):
    
    ## define a main function callable from another python script
    
    parser = argparse.ArgumentParser(description='Generate defect supercells.')
    parser.add_argument('dir_poscars',help='path to directory with unitcell POSCARs')
    parser.add_argument('json_initdef',help='json file with details to initialize defect')
    parser.add_argument('q',type=int,help='charge')
    parser.add_argument('supercell',help='supercell size')
    parser.add_argument('vacuum',type=int,help='vacuum spacing')
    parser.add_argument('--write_bulkref',help='to write bulkref?',default=False,action='store_true')
      
    ## read in the above arguments from command line
    args = parser.parse_args(args)
    
    
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
    
    
    # read in defect details from initdefect json file
    with open(args.json_initdef, 'r') as file:
        initdef = json.loads(file.read())

    ## initialize defect object
    defect = Defect(structure_bulk,structure.copy(),nvecs,args.vacuum,args.q)

    ## set the defect info (type, site, species) for each defect listed in the json file
    for d in initdef:
        initdef[d]["index_offset_n1"] = initdef[d].get("index_offset_n1",0)
        initdef[d]["index_offset_n2"] = initdef[d].get("index_offset_n2",0)
        initdef[d]["index_offset_n1n2"] = initdef[d].get("index_offset_n1n2",0)
        print (initdef[d])
        defect_site = defect.get_defect_site(initdef[d])[0]
        defect.add_defect_info(initdef[d],defect_site)  

    ## create vacancy defect(s)           
    defect.remove_atom()
    ## create substitutional defect(s)
    defect.replace_atom()
    ## create interstitial defect(s)
    defect.add_atom()     

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
    
    
if __name__ == '__main__':

    main(sys.argv[1:])    

    
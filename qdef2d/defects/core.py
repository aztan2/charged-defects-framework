import itertools
import numpy as np
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pymatgen
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.sites import PeriodicSite


class Defect():
    
    def __init__(self, structure_bulk, structure, supercell, vacuum, charge):
        
        self.structure_bulk = structure_bulk
        self.structure = structure
        self.supercell = supercell
        self.vacuum = vacuum
        self.charge = charge
        self.defect_type = []
        self.defect_site = []

   
    def parse_initdef(self, initdef):
        
        if "type" not in initdef.keys() or initdef["type"][0] not in ["v","i","s","a"]:
            raise ValueError ("unfamiliar defect type specified in initdef")
            return False
            
        elif initdef["type"][0] == "v":
            ## vacancy defect type requires the following fields:
            for field,typ in zip(["species","index"],[str,int]):
                if field not in initdef.keys():
                    raise ValueError ("missing %s entry in initdef"%field)
                    return False
                elif not isinstance(initdef[field],typ):
                    raise ValueError ("wrong entry type for %s in initdef"%field)
                    return False
            for field in ["index_offset_basis","index_offset_a","index_offset_b"]:
                if field in initdef.keys():
                    if not isinstance(initdef[field],int):
                        raise ValueError ("wrong entry type for %s in initdef"%field)
                        return False
        """                
        elif initdef["type"][0] == "i":                        
            ## interstitial defect type requires the following fields:
            for field,typ in zip(["species","index"],[str,int]):
                if field not in initdef.keys():
                    raise ValueError ("missing %s entry in initdef"%field)
                elif not isinstance(initdef[field],typ):
                    raise ValueError ("wrong entry type for %s in initdef"%field)
            for field in ["index_offset_basis","index_offset_a","index_offset_b"]:
                if field in initdef.keys():
                    if not isinstance(initdef[field],int):
                        raise ValueError ("wrong entry type for %s in initdef"%field)
            
                
        
        elif initdef["type"][0] == "i" or initdef["type"][0] == "a":
            if len(initdef["index"]) != len(initdef["species"]):
                raise ValueError ("inconsistency between index and species lists")
            for offset in ["index_offset_basis","index_offset_a","index_offset_b"]:
                if initdef[offset] == 0:
                    initdef[offset] = [0]*len(initdef["index"])
                    """
        return initdef
    
    def get_site_ind(self, site_sp, rel_ind, offset_basis_atom, offset_a, offset_b):
        
        """
        Get the absolute site index based on the relative site index 
        and any additional offsets

        Parameters
        ----------
        site_sp (str): atom species at the site of interest
        rel_ind (int): relative site index
        offset_basis_atom (int): if > 1 atom of same species in the unit cell,
                                 determine site index relative to different basis atom
        offset_a (int): # units to offset site along cell vector a
        offset_b (int): # units to offset site along cell vector b
        

        Returns
        -------
        siteind: absolute site index in supercell

        """

        ## figuring out min and max+1 indices for each species        
        elems_all = [site.specie.symbol for site in self.structure]
        elems_uniq = [el[0] for el in itertools.groupby(elems_all)]
        elem_max = np.cumsum([len(tuple(el[1])) for el in itertools.groupby(elems_all)])
        elem_min = [0] + [i for i in elem_max[:-1]]
        elem_minmax = {el: [imin,imax] for el,imin,imax in zip(elems_uniq,elem_min,elem_max)}

        ## evaluate the absolute site index
        ## I don't love this clunky method of specifying sites
        ## but I have yet to figure out a better way to do it... :(
        siteind = (rel_ind + offset_basis_atom*self.supercell[0]*self.supercell[1]
                           + offset_a*self.supercell[1] + offset_b)
        if siteind < 0:
            return(int(elem_minmax[site_sp][1]+siteind))
        else:
            return(int(elem_minmax[site_sp][0]+siteind))  
    

    def get_defect_site(self, initdef):

        siteinds = []
        
        ## this is relatively simple for vacancies and substitutionals
        ## which by definition are a pre-existing site in the structure        
        if initdef["type"][0] == "v" or initdef["type"][0] == "s":
            siteind = self.get_site_ind(initdef["species"],
                                        initdef["index"],
                                        initdef["index_offset_basis"],
                                        initdef["index_offset_a"],
                                        initdef["index_offset_b"])
            siteinds.append(siteind)
            defect_site = self.structure_bulk[siteind]

        
        ## for interstitials/adatoms, the defect site coords are defined by
        ## the center (average) of the user-provided list of sites
        elif initdef["type"][0] == "i" or initdef["type"][0] == "a":
            if len(initdef["index"]) != len(initdef["species"]):
                raise ValueError ("inconsistency between index and species lists")
            for offset in ["index_offset_basis","index_offset_a","index_offset_b"]:
                if initdef[offset] == 0:
                    initdef[offset] = [0]*len(initdef["index"])

            ## get absolute site indices
            siteind = [self.get_site_ind(sp,ind,offset_basis,offset_a,offset_b)
                           for sp,ind,offset_basis,offset_a,offset_b
                               in zip(initdef["species"],
                                      initdef["index"],
                                      initdef["index_offset_basis"],
                                      initdef["index_offset_a"],
                                      initdef["index_offset_b"])]
            siteinds.append(siteind)
            
            ## get the averaged position
            coords_ref = self.structure_bulk[siteind[0]].frac_coords
            print (coords_ref)
            defect_coords = coords_ref.copy()
            for i in siteind[1:]:
                coords = self.structure_bulk[i].frac_coords
                print (coords,coords_ref,np.array(coords[0:3])-np.array(coords_ref[0:3]))
                for j in [0,1,2]:
                    if abs(coords[j] - coords_ref[j]) > 0.501:  ## dealing with pbc 
                        if coords[j] > coords_ref[j]: coords[j] -= 1.0
                        else: coords[j] += 1.0
                print (coords)
                defect_coords += coords
            defect_coords = defect_coords/len(siteind)
            print (defect_coords)

            ## we may want to further shift the x-y position manually
            M = self.structure.lattice.matrix
            if initdef.get("shift_x") != None:
                defect_coords += float(initdef["shift_x"]) * np.linalg.inv(M)[0,:]
            if initdef.get("shift_y") != None:
                defect_coords += float(initdef["shift_y"]) * np.linalg.inv(M)[1,:]

            ## we may want to shift the z position, e.g. in the case of an adatom            
            if initdef.get("shift_z") != None:
                defect_coords += float(initdef["shift_z"]) * np.linalg.inv(M)[2,:]
            
            ## create the defect_site as a PeriodicSite object
            defect_site = PeriodicSite(initdef["species_new"],defect_coords,
                                       lattice=self.structure.lattice,
                                       coords_are_cartesian=False)
        
        return defect_site


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

    

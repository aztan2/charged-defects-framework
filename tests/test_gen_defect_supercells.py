import sys
sys.path.append('../')
import numpy as np
from pymatgen.core.operations import SymmOp
from pymatgen import Structure,Lattice
from pymatgen.analysis.structure_matcher import StructureMatcher
import unittest
import gen_defect_supercells


class Test2D(unittest.TestCase):

    def setUp(self):
        
        ## simple toy structure for preliminary testing
        
        self.a0 = 3.0
        self.c = 20.0

        self.structure = Structure(
                Lattice.from_parameters(a=self.a0, b=self.a0, c=self.c,
                                        alpha=90, beta=90, gamma=90),
                ["O","O"],
                [[0.0, 0.0, 0.1],
                 [0.5, 0.5, 0.3]])

        self.nvecs = [3,3,1]
        self.vacuum = 20
        self.q = 0
        self.structure.make_supercell(self.nvecs)
        self.structure_bulk = self.structure.copy()

        self.initdef_list = []        
        self.initdef_list.append({"def1": {"type": "vac", "species": "O", "index": 0}})
        self.initdef_list.append({"def1": {"type": "sub", "index": -1, "species": "O", 
                                           "species_new": "N"}})
        self.initdef_list.append({"def1": {"type": "sub", "index": -1, "species": "O", 
                                           "species_new": "N"},
                                  "def2": {"type": "vac", "index": -1, "species": "O", 
                                           "index_offset_n1n2": -1}})
        self.initdef_list.append({"def1": {"type": "ad", "index": [-1], "species": ["O"], 
                                           "species_new": "N", "shift_z": "3.0"}})    
        self.initdef_list.append({"def1": {"type": "int", "index": [0,0,0,0],
                                           "index_offset_n1": [0,1,0,0], 
                                           "index_offset_n2": [0,0,0,2],
                                           "index_offset_n1n2": [0,0,1,1],
                                           "species": 4*["O"], "species_new": "N"}})
        self.create_defects()

        self.siteinds_list = [[0],[17],[17,8],[[17]],[[0,3,9,11]]]  
        self.defcoords_list = [[[0.0,0.0,0.1]],[[0.833333,0.833333,0.3]],
                               [[0.833333,0.833333,0.3],[0.666667,0.666667,0.1]],
                               [[0.833333,0.833333,0.45]],[[0.166667,0.0,0.2]]]
        self.natoms_list = [17,18,17,19,19]


    def create_defects(self):
        
        ## create defects
        
        self.defects_list = []
        
        for initdef in self.initdef_list:
            
            ## initialize defect object
            defect = gen_defect_supercells.Defect(self.structure_bulk,
                                                  self.structure.copy(),
                                                  self.nvecs,self.vacuum,self.q)
            
            ## set the defect info (type, site, species) for each defect
            for d in initdef:
                initdef[d]["index_offset_n1"] = initdef[d].get("index_offset_n1",0)
                initdef[d]["index_offset_n2"] = initdef[d].get("index_offset_n2",0)
                initdef[d]["index_offset_n1n2"] = initdef[d].get("index_offset_n1n2",0)
                defect_site, siteinds = defect.get_defect_site(initdef[d])
                defect.add_defect_info(initdef[d],defect_site)  
                
            ## create defect(s)          
            defect.remove_atom()
            defect.replace_atom()
            defect.add_atom()
            self.defects_list.append(defect)
    
    
    def test_get_site_index(self):
        
        ## test get_site_index returns the correct absolute site indices

        ## initialize dummy "defect" object
        defect = gen_defect_supercells.Defect(self.structure_bulk,
                                              self.structure.copy(),
                                              self.nvecs,self.vacuum,self.q)

        for initdef,siteinds in zip(self.initdef_list,self.siteinds_list):
            for d,siteind_ref in zip(initdef,siteinds):

                if initdef[d]["type"][0] == "v" or initdef[d]["type"][0] == "s":
                    siteind = defect.get_site_ind(initdef[d]["index"],
                                                  initdef[d]["species"],
                                                  initdef[d]["index_offset_n1"],
                                                  initdef[d]["index_offset_n2"],
                                                  initdef[d]["index_offset_n1n2"])
                    self.assertEqual(siteind,siteind_ref)    
                    
                if initdef[d]["type"][0] == "a" or initdef[d]["type"][0] == "i":
                    for siteindi_ref,ind,sp,offset_n1,offset_n2,offset_n1n2 \
                        in zip(siteind_ref,
                               initdef[d]["index"],
                               initdef[d]["species"],
                               initdef[d]["index_offset_n1"],
                               initdef[d]["index_offset_n2"],
                               initdef[d]["index_offset_n1n2"]):
                        siteind = defect.get_site_ind(ind,sp,offset_n1,
                                                      offset_n2,offset_n1n2)
                        self.assertEqual(siteind,siteindi_ref)


    def test_defcoords(self):
        
        ## test that the defect(s) are created at the correct position
        ## essentially tests the get_defect_site function

        for defects,defcoords in zip(self.defects_list,self.defcoords_list):
            for defect_site,defcoord in zip(defects.defect_site,defcoords):
                for j in range(3):
                    self.assertAlmostEqual(defect_site.frac_coords[j]%1%1,
                                           defcoord[j]%1%1,places=6)  


    def test_natoms(self):
        
        ## test that the defect creation results in the correct number of atoms
        
        for defect,natoms in zip(self.defects_list,self.natoms_list):
            self.assertEqual(defect.structure.num_sites,natoms)
        

class Test2D_WSe2(Test2D):

    def setUp(self):
        
        ## hexagonal WSe2 unitcell
        
        self.a0 = 3.287596
        self.c = 23.360843

        self.structure = Structure(
                Lattice.from_parameters(a=self.a0, b=self.a0, c=self.c,
                                        alpha=90, beta=90, gamma=120),
                ["W","Se","Se"],
                [[0.0, 0.0, 0.5],
                 [0.333333, 0.666667, 0.571091],
                 [0.333333, 0.666667, 0.428909]])

        self.nvecs = [4,4,1]
        self.vacuum = 20
        self.q = 0
        self.structure.make_supercell(self.nvecs)
        self.structure_bulk = self.structure.copy()

        self.initdef_list = []        
        self.initdef_list.append({"def1": {"type": "vac", "species": "Se", "index": 0}})
        self.initdef_list.append({"def1": {"type": "sub", "index": -1, "species": "W", 
                                           "species_new": "Re"}})
        self.initdef_list.append({"def1": {"type": "sub", "index": -1, "species": "W", 
                                           "species_new": "Re"},
                                  "def2": {"type": "vac", "index": -1, "species": "Se", 
                                           "index_offset_n1n2": -1}})
        self.initdef_list.append({"def1": {"type": "ad-W", "index": [-1], "species": ["W"], 
                                           "species_new": "Re", "shift_z": "3.36"}})    
        self.initdef_list.append({"def1": {"type": "int-hex", "index": [0,1,0],
                                           "index_offset_n1": [1,1,0], 
                                           "species": 3*["W"], "species_new": "Re"}})
        self.create_defects()

        self.siteinds_list = [[16],[15],[15,31],[[15]],[[4,5,0]]]        
        self.defcoords_list = [[[0.083333,0.166667,0.571091]],[[0.75,0.75,0.5]],
                               [[0.75,0.75,0.5],[0.833333,0.916667,0.571091]],
                               [[0.75,0.75,0.643830]],[[0.166667,0.083333,0.5]]]
        self.natoms_list = [47,48,47,49,49]              
        
    
if __name__ == '__main__':

  
    suite = unittest.TestLoader().loadTestsFromTestCase(Test2D)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(Test2D_WSe2)
    unittest.TextTestRunner(verbosity=2).run(suite)
    

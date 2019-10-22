import sys
sys.path.append('../')
import numpy as np
from pymatgen.core.operations import SymmOp
from pymatgen import Structure,Lattice
from pymatgen.analysis.structure_matcher import StructureMatcher
import unittest
import gen_unitcell_2d


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
        
        self.structure_rot = Structure(
                Lattice.from_parameters(a=self.a0, b=self.c, c=self.a0,
                                        alpha=90, beta=90, gamma=90),
                ["O","O"],
                [[0.0, 0.1, 0.0],
                 [0.5, 0.3, 0.5]])
        self.zaxis = 'b'
        
        self.structure_bulk = Structure(
                Lattice.from_parameters(a=self.a0, b=self.a0, c=self.c,
                                        alpha=90, beta=90, gamma=90),
                ["O","O","O","O"],
                [[0.0, 0.0, 0.15],
                 [0.5, 0.5, 0.35],
                 [0.0, 0.0, 0.65],
                 [0.5, 0.5, 0.85]])
        self.slabmin = 0.0
        self.slabmax = 0.5


    def test_align_axis_matcher(self):
        
        ## test align_axis function using pymatgen's StructureMatcher

        structure_rot2 = gen_unitcell_2d.align_axis(self.structure_rot, self.zaxis)

        matcher = StructureMatcher()
        self.assertTrue(matcher.fit(self.structure, structure_rot2))


    def test_align_axis_lattice(self):
        
        ## test align_axis function by comparing lattice vectors
        
        structure_rot2 = gen_unitcell_2d.align_axis(self.structure_rot, self.zaxis)
        
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(self.structure.lattice._matrix[i][j],
                                       structure_rot2.lattice._matrix[i][j],places=8)


    def test_align_axis_coords(self):
        
        ## test align_axis function by comparing atomic coords
        
        structure_rot2 = gen_unitcell_2d.align_axis(self.structure_rot, self.zaxis)
        
        for s1,s2 in zip(self.structure.sites,structure_rot2.sites):
            for j in range(3):
                self.assertAlmostEqual(s1._fcoords[j]%1%1,s2._fcoords[j]%1%1,places=8)        
        

    def test_center_slab(self):
        
        ## test center_slab function
        
        for shift in [0.0,0.1,0.5,-0.2,2.0]:
            
            struct = self.structure.copy()
            struct.translate_sites(range(self.structure.num_sites), (0, 0, shift))
            struct = gen_unitcell_2d.center_slab(struct)
            slab_center = np.average([s._fcoords[2] for s in struct.sites])  
            
            self.assertAlmostEqual(slab_center,0.5,places=8)


    def test_layer_from_bulk(self):
        
        ## test layer_from_bulk function
        
        struct_layer = gen_unitcell_2d.layer_from_bulk(self.structure_bulk,
                                                       self.slabmin,self.slabmax)

        for s1,s2 in zip(gen_unitcell_2d.center_slab(self.structure).sites,
                         gen_unitcell_2d.center_slab(struct_layer).sites):
            for j in range(3):
                self.assertAlmostEqual(s1._fcoords[j]%1%1,s2._fcoords[j]%1%1,places=8)  
                

class Test2D_WSe2(Test2D):

    def setUp(self):
        
        ## hexagonal WSe2 unitcell
        
        self.a0 = 3.325612
        self.c = 17.527085

        self.structure = Structure(
                Lattice.from_parameters(a=self.a0, b=self.a0, c=self.c,
                                        alpha=90, beta=90, gamma=120),
                ["W","Se","Se"],
                [[0.0, 0.0, 0.0],
                 [0.333333, 0.666667, 0.095876],
                 [0.333333, 0.666667, 0.904124]])
       
        self.structure_rot = Structure(
                Lattice.from_parameters(a=self.c, b=self.a0, c=self.a0,
                                        alpha=60, beta=90, gamma=90),
                ["W","Se","Se"],
                [[0.0, 0.0, 0.0],
                 [0.095876, 0.666667, -0.333333],
                 [0.904124, 0.666667, -0.333333]])    
        self.zaxis = 'a'

        self.structure_bulk = Structure(
                Lattice.from_parameters(a=self.a0, b=self.a0, c=self.c,
                                        alpha=90, beta=90, gamma=120),
                ["W","W","Se","Se","Se","Se"],
                [[0.0, 0.0, 0.75],
                 [0.333333, 0.666667, 0.25],
                 [0.0, 0.0, 0.345876],
                 [0.333333, 0.666667, 0.845876],
                 [0.333333, 0.666667, 0.654124],
                 [0.0, 0.0, 0.154124]])
        self.slabmin = 0.5
        self.slabmax = 1.0        


class Test2D_SnS(Test2D):

    def setUp(self):
        
        ## SnS unitcell
        
        self.a0 = 4.442511
        self.a1 = 4.023972
        self.c = 11.432652

        self.structure = Structure(
                Lattice.from_parameters(a=self.a0, b=self.a1, c=self.c,
                                        alpha=90, beta=90, gamma=90),
                ["Sn","Sn","S","S"],
                [[0.873976, 0.250000, 0.121238],
                 [0.373976, 0.750000, 0.378762],
                 [0.480133, 0.750000, 0.150170],
                 [0.980133, 0.250000, 0.349830]])

        self.structure_rot = Structure(
                Lattice.from_parameters(a=self.c, b=self.a1, c=self.a0,
                                        alpha=90, beta=90, gamma=90),
                ["Sn","Sn","S","S"],
                [[0.121238, 0.250000, 0.126024],
                 [0.378762, 0.750000, 0.626024],
                 [0.150170, 0.750000, 0.519867],
                 [0.349830, 0.250000, 0.019867]])
        self.zaxis = 'a'
        
        self.structure_bulk = Structure(
                Lattice.from_parameters(a=self.a0, b=self.a1, c=self.c,
                                        alpha=90, beta=90, gamma=90),
                ["Sn","Sn","Sn","Sn","S","S","S","S"],
                [[0.873976, 0.250000, 0.121238],
                 [0.126024, 0.750000, 0.878762],
                 [0.373976, 0.750000, 0.378762],
                 [0.626024, 0.250000, 0.621238],
                 [0.480133, 0.750000, 0.150170],
                 [0.519867, 0.250000, 0.849830],
                 [0.980133, 0.250000, 0.349830],
                 [0.019867, 0.750000, 0.650170]])
        self.slabmin = 0.0
        self.slabmax = 0.5 
        
    
if __name__ == '__main__':

  
    suite = unittest.TestLoader().loadTestsFromTestCase(Test2D)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(Test2D_WSe2)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    suite = unittest.TestLoader().loadTestsFromTestCase(Test2D_SnS)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    
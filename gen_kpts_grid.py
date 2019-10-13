import sys
import os
import math
import argparse
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pymatgen
from pymatgen.io.vasp.inputs import Poscar,Kpoints


def automatic_density_2d(structure, kppa, force_gamma=False):

    comment = "automatically generated KPOINTS with 2d grid density = " + \
        "%.0f per reciprocal atom" % kppa
#    if math.fabs((math.floor(kppa ** (1 / 3) + 0.5)) ** 3 - kppa) < 1:
    if math.fabs((math.floor(kppa ** (1 / 2) + 0.5)) ** 2 - kppa) < 1:
        kppa += kppa * 0.01
    latt = structure.lattice
    lengths = latt.abc
    ngrid = kppa / structure.num_sites
#    mult = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1 / 3)
    mult = (ngrid * lengths[0] * lengths[1]) ** (1 / 2)

    num_div = [int(math.floor(max(mult / l, 1))) for l in lengths]
    num_div[2] = 1  ## force only 1 kpt in c direction

    is_hexagonal = latt.is_hexagonal()

    has_odd = any([i % 2 == 1 for i in num_div])
    if has_odd or is_hexagonal or force_gamma:
        style = Kpoints.supported_modes.Gamma
    else:
        style = Kpoints.supported_modes.Monkhorst

    return Kpoints(comment, 0, style, [num_div], [0, 0, 0])


def main(args):
    
    parser = argparse.ArgumentParser(description='Generate KPOINTS')
    parser.add_argument('--kppa',type=int,help='kpt density (pra)',default=440)
      
    ## read in the above arguments from command line
    args = parser.parse_args(args)
    
    ## the bash script already put us in the appropriate subdirectory
    dir_sub = os.getcwd()
    
    poscar = Poscar.from_file(os.path.join(dir_sub,"POSCAR"))
    kpts = automatic_density_2d(poscar.structure, args.kppa, force_gamma=False) 
    
    Kpoints.write_file(kpts, os.path.join(dir_sub,"KPOINTS"))


if __name__ == '__main__':
    
    main(sys.argv[1:]) 
  
    
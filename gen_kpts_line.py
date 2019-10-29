import sys
import os
import argparse
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pymatgen
from pymatgen.io.vasp.inputs import Poscar, Kpoints
from pymatgen.symmetry.bandstructure import HighSymmKpath


def remove_z_kpoints(kpts_line):
    
    """
    Strips all k-points from the standard linemode k-points object
    that include a non-zero z-component,
    since these are not relevant for 2D materials and slabs.
    
    Parameters
    ----------
    kpts_line: Kpoint object describing the high-symmetry k-path

    Returns
    -------
    kpts_2d: Kpoint object describing the high-symmetry k-path for 2D structure
    
    """
    
    num_branches = int(len(kpts_line.kpts)/2)

    kpts_2d_ind = []    
    for i in range(num_branches):
        ## alternating kpts define the start and end points for each branch
        ## only keep those branches in which both endpoints have zero z-component
        if kpts_line.kpts[2*i][2] == 0.0 and kpts_line.kpts[2*i+1][2] == 0.0:
            kpts_2d_ind.append(2*i) 
            kpts_2d_ind.append(2*i+1)  

    kpts_2d = kpts_line
    kpts_2d.kpts = [kpts_line.kpts[i] for i in kpts_2d_ind]
    kpts_2d.labels = [kpts_line.labels[i] for i in kpts_2d_ind]  
    
    return kpts_2d


def get_ibzkpts(ibzkpt_path="../"):
    
    """
    Gets the uniform kpt mesh from IBZKPT file.

    Parameters
    ----------
    ibzkpt_path: path to IBZKPT file from a previous structural relaxation run
                 (default assumes one directory up)

    Returns
    -------
    ibz_lines : lines of IBZKPT file
    n_ibz_kpts : number of symmetrically unique IBZ kpoints

    """
        
    ibz_lines = open(os.path.join(ibzkpt_path,"IBZKPT")).readlines()
    for i, line in enumerate(ibz_lines):
        if "Tetrahedra" in line:
            ibz_lines = ibz_lines[:i]
            break
    n_ibz_kpts = int(ibz_lines[1].split()[0])
    
    return ibz_lines, n_ibz_kpts


def get_kpts_line_explicit(kpts_line, ndiv):
    
    """
    Explicitly lists all points along a given high-symmetry k-path

    Parameters
    ----------
    kpts_line: Kpoint object describing the high-symmetry k-path
    ndiv : (int) number of divisions along high-symmetry lines

    Returns
    -------
    abs_path : list of all points along the given high-symmetry k-path
        
    """         
    
    num_branches = int(len(kpts_line.kpts)/2)
    
    abs_path = []
    for i in range(num_branches):
        ## alternating kpts define the start and end points for each branch
        start_kpt = [kpts_line.kpts[2*i][0],kpts_line.kpts[2*i][1],
                     kpts_line.kpts[2*i][2],str(kpts_line.labels[2*i])]
        end_kpt = [kpts_line.kpts[2*i+1][0],kpts_line.kpts[2*i+1][1],
                   kpts_line.kpts[2*i+1][2],str(kpts_line.labels[2*i+1])]
#        print (start_kpt,end_kpt)
        
        increments = [(end_kpt[0] - start_kpt[0]) / ndiv,
                      (end_kpt[1] - start_kpt[1]) / ndiv,
                      (end_kpt[2] - start_kpt[2]) / ndiv]
        
        abs_path.append([str(k) for k in start_kpt[:3]] + ['0', start_kpt[-1]])
        ## explicitly write out all the intermediate kpts along the branch
        for n in range(1, ndiv):
            abs_path.append([str(start_kpt[0] + increments[0] * n),
                             str(start_kpt[1] + increments[1] * n),
                             str(start_kpt[2] + increments[2] * n), '0'])
        abs_path.append([str(k) for k in end_kpt[:3]] + ['0', end_kpt[-1]])            

    return abs_path
            

def main(args):
    
    """
    Uses pymatgen symmetry functions to determine the high-symmetry k-path.
    If 2D structure specified, removes k-points with non-zero z-component.
    
    Parameters
    ----------
    structure (Structure): structure for determining k-path
    ndiv (int): number of divisions along high-symmetry lines
    dim (int): 2 for a 2D material, 3 for a 3D material.
    
    Returns
    -------
    kpts_line : Kpoint object
        
    """
   
    parser = argparse.ArgumentParser(description='Generate linemode KPOINTS')
    parser.add_argument('--ndiv',type=int,default=20,
                        help='number of divisions along high-symmetry lines')
    parser.add_argument('--dim',type=int,default=2,
                        help='dimensionality of structure')
    parser.add_argument('--for_scan_hse',default=False,action='store_true',
                        help='generate special KPOINTS for SCAN/HSE?')
      
    ## read in the above arguments from command line
    args = parser.parse_args(args)
    
    ## the bash script already put us in the appropriate subdirectory
    dir_sub = os.getcwd()
    
    poscar = Poscar.from_file(os.path.join(dir_sub,"POSCAR"),
                              check_for_POTCAR=False, read_velocities=False)
    structure = poscar.structure

    ## use pymatgen symmetry functions to determine the high-symmetry k-path
    kpts_line = Kpoints.automatic_linemode(args.ndiv, HighSymmKpath(structure))
    ## if 2D structure specified, remove k-points with non-zero z-component
    if args.dim == 2:
        kpts_line = remove_z_kpoints(kpts_line)  
    
    if not args.for_scan_hse:
        ## write out regular KPOINT file for bandstructure calc. with PBE
        kpts_line.write_file(os.path.join(dir_sub,"KPOINTS_bands"))
        
    else:
        ## get regular grid kpoints from IBZKPT file
        ibz_lines, n_ibz_kpts = get_ibzkpts()
        ## get explicit kpoints between each high-symmetry point
        abs_path = get_kpts_line_explicit(kpts_line, args.ndiv)
        ## write out special KPOINT file for bandstructure calc. with SCAN/HSE
        with open(os.path.join(dir_sub,"KPOINTS_bands_SCAN"),'w') as f:
            f.write('Special KPOINTS file\n')
            f.write('{}\n'.format(n_ibz_kpts + len(abs_path)))
            f.write('Reciprocal Lattice\n')
            for line in ibz_lines[3:]:
                f.write(line)
            for point in abs_path:
                f.write('{}\n'.format(' '.join(point)))


if __name__ == '__main__':

    """
    Writes a KPOINTS file for band structure calculations. Does
    not use the typical linemode syntax for NSCF calculations,
    but uses the IBZKPT + high-symmetry path syntax described in
    http://cms.mpi.univie.ac.at/wiki/index.php/Si_bandstructure
    so that SCF calculations can be performed.
    
    """
    
    main(sys.argv[1:]) 

            
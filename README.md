# Framework for first-principles calculations of charged dopants and defects in 2D materials

I am developing a Python-based workflow for efficient high-throughput calculations of charged point defects in 2D materials to create a new database of defects in 2D and layered bulk materials. The development of this workflow will enable the rapid and robust evaluation of defect structures, formation energies, and electronic properties for different types of defects in a wide range of 2D materials.

This repository contains a collection of python and bash scripts for setting up and post-processing of charged defects calculations in 2D materials. The scripts are set up to work with VASP input/output files and the University of Florida's high performance computing cluster, HiPerGator. 

This is currently a work in progress by Anne Marie Tan. 

## Python scripts: 

Note: These python scripts make use of the following python packages: `numpy`, `pymatgen`, `pandas`, `matplotlib`, `argparse` which need to be installed beforehand.

* `gen_defect_supercells.py`: generates a defect supercell given supercell size and vacuum spacing and write out the `POSCAR` and `defectproperty.json` files. ** the type and position of the defect to create is still hardcoded and requires the user to make changes in the script itself **
* `gen_incar.py`: generates an `INCAR` file with some standard settings depending on the user-specified runtype (relaxation/dos), functional (PBE/SCAN+rVV10).
* `gen_kpts_grid.py`: generates a `KPOINTS` file with standard Monkorst-Pack grid and user-specified kpoint density.
* `gen_submit.py`: generates a SLURM submison script for submitting jobs on HiPerGator. User can specify the queue, nodes, memory, time limit to be requested.

* `gen_SPHInX_input_file.py`: generates the input file required for the [Freysoldt 2D charge correction scheme](https://doi.org/10.1103/PhysRevB.97.205425).
* `get_alignment_correction_2d.py`: iteratively applies the [sxdefectalign2d script](https://sxrepo.mpie.de/projects/sphinx-add-ons/files) until the optimal charge position and correction energy are obtained.
* `get_alignment_correction.py`: post-processing of the [sxdefectalign script](https://sxrepo.mpie.de/projects/sphinx-add-ons/files) to determine an averaged alignment correction for the 3D case. (this script is old and has not been tested recently)

* `parse_energies.py`: parses the total energies from VASP `vasprun.xml` files, stores in a pandas dataframe, and writes out to an excel file.
* `parse_corrections.py`: parses the correction energies computed by sxdefectalign2d, adds to the pandas dataframe, and writes out to an excel file.


## Bash scripts: 

* `setup_defects.sh`: creates new sub-directories in the charge/supercell/vacuum pattern for a given defect, and runs `gen_defect_supercells.py`, `gen_incar.py`, `gen_kpts_grid.py`, and `gen_submit.py` to generate all the required input files for a VASP calculation. This script assumes that the parent directory from which this is run contains the appropriate unitcell `POSCAR` and `POTCAR` files.
* `apply_correction_loop.sh`: runs `gen_SPHInX_input_file.py` followed by `get_alignment_correction_2d.py` in all relevant sub-directories (ignoring neutral defects)
* `restart_SCAN.sh`: creates new sub-directories in the charge/supercell/vacuum pattern and prepares the input files to run VASP calculations with SCAN+rVV10 functional, starting from PBE-relaxed `CONTCAR` (and `WAVECAR` if available).
* `restart_soc.sh`: creates new sub-directories and prepares the input files to run a single-point density of states VASP calculation with spin-orbit coupling, starting from the appropriate relaxed `CONTCAR`.


## Authors:
Anne Marie Z. Tan


## How to cite:
BibTex entry for this Github repository::

```
   @misc{charged-defects-framework,
     title        = {Framework for first-principles calculations of charged dopants and defects in 2D materials},
     author       = {A. M. Z Tan},
     year         = 2019,
     publisher    = {GitHub},
     journal      = {GitHub repository},
     howpublished = {\url{https://github.com/aztan2/charged-defects-framework}},
     url          = {https://github.com/aztan2/charged-defects-framework},
     doi          = {}
   }
```



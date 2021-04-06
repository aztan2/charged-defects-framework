# Framework for first-principles calculations of charged dopants and defects in 2D materials

**qdef2d** [tentative name, work-in-progress] is a python package that aims to help automate and accelerate the setting up and post-processing of charged defect calculations in 2D materials. This will enable the rapid and robust evaluation of defect structures, formation energies, and electronic properties for different types of defects in a wide range of 2D materials to create a new database of defects in 2D and layered bulk materials. The modules and scripts are currently set up to work with VASP input/output files and the University of Florida's high performance computing cluster, HiPerGator. 


## Requirements:
* python3 environment with the following packages installed: `numpy`, `pymatgen`, `pandas` (+ `openpyxl`), `matplotlib` 
* download the [sxdefectalign2d script](https://sxrepo.mpie.de/projects/sphinx-add-ons/files) which we will be using to apply the [Freysoldt-Neugebauer 2D charge correction scheme](https://doi.org/10.1103/PhysRevB.97.205425).


## qdef2d package:

```
qdef2d
|-- io
|   |-- vasp
|   |   |-- incar.py
|   |   |-- kpoints.py
|   |   |-- submit.py
|   |   |-- parse_energies.py
|   |-- database
|       |-- database_entry.py
|-- defects
|   |-- core.py
|   |-- gen_defect_supercell.py
|   |-- setup_defect_calcs.py
|   |-- calc_Eform_uncorr.py
|   |-- calc_Eform_corr.py
|   |-- corrections
|       |-- SPHInX_input_file.py
|       |-- alignment_correction_2d.py
|       |-- apply_corrections_2d.py
|       |-- parse_corrections.py
|-- logging.py
|-- osutils.py
|-- slabutils.py
```

* `incar.py`, `kpoints.py`, `submit.py`: functions for generating VASP input files, job submission script
* `parse_energies.py`: function/script to parse total energies from VASP OUTCARs and save into a pandas dataframe.
* `database_entry.py`: functionalities related to creating, manipulating, and reading simple database entries [should be replaced with interface to actual mongodb database, e.g. on MaterialsWeb]
* `core.py`: defines Defect object, includes functionalities for creating different types of defects
* `gen_defect_supercell.py`, `setup_defect_calcs.py`: functions/scripts to generate defect supercell, VASP input files for defect supercell calculations
* `calc_Eform_uncorr.py`, `calc_Eform_corr.py`: functions/scripts to evaluate the un-corrected and corrected defect formation energies and save into a pandas dataframe.
* `SPHInX_input_file.py`, `alignment_correction_2d.py`, `apply_corrections_2d.py`, `parse_corrections.py`: functions/scripts to generate the input file for the charge correction, apply the charge correction, and parse the results into a pandas dataframe.
* `logging.py`, `osutils.py`, `slabutils.py`: additional utility functions

Additional documentation can be found in the source codes, accessed via the `help()` function, or with the `--h` flag from the command line.


## Jupyter tutorials:

These jupyter notebooks demonstrate how to make use of the modules and scripts in **qdef2d** to set up, execute, and analyze the results of DFT calculations to obtain various properties of 2D materials. They are meant to be executed on HiPerGator, where pre-calculated examples can also be found. 

* `tutorial_2D_pristine.ipynb`: general tutorial about how to calculate some basic properties (lattice constant, band structure, dielectric tensor) of pristine 2D materials.
* `tutorial_2D_defects.ipynb`: tutorial about how to correctly evaluate the formation energies of charged defects in 2D materials.


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

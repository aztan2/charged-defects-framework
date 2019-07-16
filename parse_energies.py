import os
import argparse
import logging
import time
import numpy as np
import pandas as pd
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun


def listdironly(path):
    
    return [d for d in os.listdir(path) if os.path.isdir(path+d) == True]


def joinpath(path,*paths):
    
    ## force the paths to have forward slash (unix-style)
    return os.path.join(path,*paths).replace('\\','/')


if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description='Parse total energies from OUTCARs.')
    parser.add_argument('path',help='path to the directory containing all the output files')
    parser.add_argument('path_ref',help='path to the directory containing all the reference output files')
    parser.add_argument('xlfile',help='excel filename to save the dataframe to')
    parser.add_argument('--logfile',help='logfile to save output to')
       
    ## read in the above arguments from command line
    args = parser.parse_args()
    
    ## set up logging
    if args.logfile:
        logging.basicConfig(filename=args.logfile,filemode='w',
                            format='%(levelname)s:%(message)s',level=logging.DEBUG)    
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(logging.Formatter('%(levelname)s:%(message)s'))
        logging.getLogger('').addHandler(console)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.DEBUG)

    qs = listdironly(args.path)
    
    writer = pd.ExcelWriter(joinpath(args.path,args.xlfile))

    time0 = time.time()    

    ## set up dataframe for neutral defect first
    if 'charge_0' not in qs:
        logging.warning("can't find output files for neutral defect")
    else:
        df0 = pd.DataFrame(columns = ['vacuum',
                                      'supercell',
                                      'N',
                                      '1/N',
                                      'E_def',
                                      'E_bulk'])
    
        for cell in listdironly(joinpath(args.path,'charge_0','')):
            for vac in listdironly(joinpath(args.path,'charge_0',cell,'')):   
                logging.info("parsing neutral %s %s"%(cell,vac))
                
                folder = joinpath(args.path,'charge_0',cell,vac,'')
#                folder = joinpath(args.path,'charge_0',cell,vac,'dos','')
                folder_ref = joinpath(args.path_ref,'charge_0',cell,vac,'bulkref','')
#                folder_ref = joinpath(args.path_ref,'charge_0',cell,vac,'bulkref','dos','')
                vr_file = joinpath(folder,'vasprun.xml')
                vr_ref_file = joinpath(folder_ref,'vasprun.xml')
                
                if not os.path.exists(vr_file):
                    logging.warning("%s file does not exist"%vr_file)
                    
                elif not os.path.exists(vr_ref_file):
                    logging.warning("%s file does not exist"%vr_ref_file)
                    
                else:
                    natoms = np.sum(Poscar.from_file(joinpath(folder_ref,'POSCAR')).natoms)
                    vr = Vasprun(vr_file)
                    vr_ref = Vasprun(vr_ref_file)
                    
                    if not vr.converged:
                        logging.warning("Vasp calculation in %s may not be converged"
                                        %folder)                    
                    if not vr_ref.converged:
                        logging.warning("Vasp calculation in %s may not be converged"
                                        %(joinpath(folder,'bulkref','')))
                    
                    df0.loc[len(df0)] = [vac,
                                         cell,
                                         natoms,
                                         1/natoms,
                                         vr.final_energy,
                                         vr_ref.final_energy]
                    ## add the parsing of chem pot, vbm
                    
        df0.sort_values(['vacuum','N'],inplace=True)
        df0.to_excel(writer,'charge_0')
    

    ## modify dataframe for charged defects
    for q in [qi for qi in qs if qi != 'charge_0']:
        df = df0.copy(deep=True)
        
        for cell in listdironly(joinpath(args.path,q,'')):
            for vac in listdironly(joinpath(args.path,q,cell,'')):
                logging.info("parsing %s %s %s"%(q,cell,vac))
                
                folder = joinpath(args.path,q,cell,vac,'')
#                folder = joinpath(args.path,q,cell,vac,'dos','')
                vr_file = joinpath(folder,'vasprun.xml')
                
                if not os.path.exists(vr_file):
                    logging.warning("%s file does not exist"%vr_file)
                    
                else:
                    vr = Vasprun(vr_file)                    
                    if not vr.converged:
                        logging.warning("Vasp calculation in %s may not be converged"
                                        %folder)                       
                    df.loc[(df['vacuum'] == vac) & 
                           (df['supercell'] == cell),'E_def'] = vr.final_energy
        
        df.to_excel(writer, q)

    writer.save()
    
    logging.debug("Total time taken (s): %.2f"%(time.time()-time0))
            
    
    

import os
import time
import argparse
import numpy as np
import pandas as pd
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from qdef2d import osutils, logging


def parse(path_def,path_ref,xlfile,soc=False,logfile=None):

    
    """ 
    Parse total energies from OUTCARs and save into a pandas dataframe.
    
    Parameters
    ----------
    path_def (str): path to the directory containing all the defect output files
    path_ref (str): path to the directory containing all the reference output files
    xlfile (str): excel filename to save the dataframe to
    [optional] soc (bool): whether or not to look in soc(dos) subdirectory. Default=False.
    [optional] logfile (str): logfile to save output to
    
    """
    
    ## set up logging
    if logfile:
        myLogger = logging.setup_logging(logfile)
    else:
        myLogger = logging.setup_logging()
        

    qs = osutils.listdironly(path_def)
    
    writer = pd.ExcelWriter(os.path.join(path_def,xlfile))

    time0 = time.time()    

    ## set up dataframe for neutral defect first
    if 'charge_0' not in qs:
        myLogger.warning("can't find output files for neutral defect")
    else:
        df0 = pd.DataFrame(columns = ['vacuum',
                                      'supercell',
                                      'N',
                                      '1/N',
                                      'E_def',
                                      'E_bulk'])
    
        for cell in osutils.listdironly(os.path.join(path_def,'charge_0')):
            for vac in osutils.listdironly(os.path.join(path_def,'charge_0',cell)):   
                myLogger.info("parsing neutral %s %s"%(cell,vac))

                subdir_def = os.path.join(path_def,'charge_0',cell,vac)
                subdir_ref = os.path.join(path_ref,'charge_0',cell,vac)

                if soc:
                    dir_soc = [dirname for dirname in osutils.listdironly(subdir_def) if "soc" in dirname]
                    if len(dir_soc) == 0:
                        myLogger.info("cannot find a soc subdirectory")
                    if len(dir_soc) > 1:
                        myLogger.info("multiple possible soc subdirectories found")
                    if len(dir_soc) == 1:
                        subdir_def = os.path.join(subdir_def,dir_soc[0])
                        subdir_ref = os.path.join(subdir_ref,dir_soc[0])
                        myLogger.info("parsing soc subdirectory")
                    
                if os.path.exists(subdir_def) and 'restart' in osutils.listdironly(subdir_def):
                    subdir_def = os.path.join(subdir_def,'restart')
                    myLogger.info("parsing restart subdirectory")
                                        
                vr_file = os.path.join(subdir_def,'vasprun.xml')
                vr_ref_file = os.path.join(subdir_ref,'vasprun.xml')
                
                if not os.path.exists(vr_file):
                    myLogger.warning("%s file does not exist"%vr_file)
                    
                elif not os.path.exists(vr_ref_file):
                    myLogger.warning("%s file does not exist"%vr_ref_file)
                    
                else:
                    natoms = np.sum(Poscar.from_file(os.path.join(subdir_ref,'POSCAR')).natoms)
                    vr = Vasprun(vr_file)
                    vr_ref = Vasprun(vr_ref_file)
                    
                    if not vr.converged:
                        myLogger.warning("VASP calculation in %s may not be converged"%subdir_def)                    
                    
                    df0.loc[len(df0)] = [vac,
                                         cell,
                                         natoms,
                                         1/natoms,
                                         vr.final_energy,
                                         vr_ref.final_energy]
                    
                    
        df0.sort_values(['vacuum','N'],inplace=True)
        df0.to_excel(writer,'charge_0', index=False)
    

    ## modify dataframe for charged defects
    for q in [qi for qi in qs if qi != 'charge_0']:
        df = df0.copy(deep=True)

        for cell in osutils.listdironly(os.path.join(path_def,q)):
            for vac in osutils.listdironly(os.path.join(path_def,q,cell)):   
                myLogger.info("parsing %s %s %s"%(q,cell,vac))

                subdir_def = os.path.join(path_def,q,cell,vac)

                if soc:
                    dir_soc = [dirname for dirname in osutils.listdironly(subdir_def) if "soc" in dirname]
                    if len(dir_soc) == 0:
                        myLogger.info("cannot find a soc subdirectory")
                    if len(dir_soc) > 1:
                        myLogger.info("multiple possible soc subdirectories found")
                    if len(dir_soc) == 1:
                        subdir_def = os.path.join(subdir_def,dir_soc[0])
                        myLogger.info("parsing soc subdirectory")
                    
                if os.path.exists(subdir_def) and 'restart' in osutils.listdironly(subdir_def):
                    subdir_def = os.path.join(subdir_def,'restart')
                    myLogger.info("parsing restart subdirectory")
                                        
                vr_file = os.path.join(subdir_def,'vasprun.xml')
                
                if not os.path.exists(vr_file):
                    myLogger.warning("%s file does not exist"%vr_file)
                    
                else:
                    vr = Vasprun(vr_file)
                    
                    if not vr.converged:
                        myLogger.warning("VASP calculation in %s may not be converged"%subdir_def) 
                        
                    df.loc[(df['vacuum'] == vac) & 
                           (df['supercell'] == cell),'E_def'] = vr.final_energy
        
        df.to_excel(writer, q, index=False)

    writer.save()
    
    myLogger.debug("Total time taken (s): %.2f"%(time.time()-time0))
    
    
if __name__ == '__main__':
    

    ## this script can also be run directly from the command line
    parser = argparse.ArgumentParser(description='Parse total energies from OUTCARs.')
    parser.add_argument('path_def',help='path to the directory containing all the defect output files')
    parser.add_argument('path_ref',help='path to the directory containing all the reference output files')
    parser.add_argument('xlfile',help='excel filename to save the dataframe to')
    parser.add_argument('--soc',help='whether or not to look in soc(dos) subdirectory',default=False,action='store_true')
    parser.add_argument('--logfile',help='logfile to save output to')

    ## read in the above arguments from command line
    args = parser.parse_args()
    
    parse(args.path_def, args.path_ref, args.xlfile, args.soc, args.logfile) 
             

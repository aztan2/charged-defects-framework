import os
import time
import argparse
import numpy as np
import pandas as pd
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
import myutils


def main(args):
    
    ## define a main function callable from another python script

    parser = argparse.ArgumentParser(description='Parse total energies from OUTCARs.')
    parser.add_argument('path',help='path to the directory containing all the output files')
    parser.add_argument('path_ref',help='path to the directory containing all the reference output files')
    parser.add_argument('xlfile',help='excel filename to save the dataframe to')
    parser.add_argument('--soc',help='whether or not to look in soc(dos) subdirectory',default=False,action='store_true')
    parser.add_argument('--logfile',help='logfile to save output to')
       
    ## read in the above arguments from command line
    args = parser.parse_args(args)
    
    ## set up logging
    if args.logfile:
        myLogger = myutils.setup_logging(args.logfile)
    else:
        myLogger = myutils.setup_logging()
        

    qs = myutils.listdironly(args.path)
    
    writer = pd.ExcelWriter(myutils.joinpath(args.path,args.xlfile))

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
    
        for cell in myutils.listdironly(myutils.joinpath(args.path,'charge_0','')):
            for vac in myutils.listdironly(myutils.joinpath(args.path,'charge_0',cell,'')):   
                myLogger.info("parsing neutral %s %s"%(cell,vac))

                folder = myutils.joinpath(args.path,'charge_0',cell,vac,'')
#                folder_ref = myutils.joinpath(args.path_ref,'charge_0',cell,vac,'bulkref','')
                folder_ref = myutils.joinpath(args.path_ref,'charge_0',cell,vac,'')

                if args.soc:  
                    folder = myutils.joinpath(folder,'dos','')
                    folder_ref = myutils.joinpath(folder_ref,'dos','')
                    myLogger.info("parsing dos subdirectory")
                if os.path.exists(folder) and 'restart' in myutils.listdironly(folder):
                    folder = myutils.joinpath(folder,'restart','')
                    myLogger.info("parsing restart subdirectory")
                                        
                vr_file = myutils.joinpath(folder,'vasprun.xml')
                vr_ref_file = myutils.joinpath(folder_ref,'vasprun.xml')
                
                if not os.path.exists(vr_file):
                    myLogger.warning("%s file does not exist"%vr_file)
                    
                elif not os.path.exists(vr_ref_file):
                    myLogger.warning("%s file does not exist"%vr_ref_file)
                    
                else:
                    natoms = np.sum(Poscar.from_file(myutils.joinpath(folder_ref,'POSCAR')).natoms)
                    vr = Vasprun(vr_file)
                    vr_ref = Vasprun(vr_ref_file)
                    
                    if not vr.converged:
                        myLogger.warning("Vasp calculation in %s may not be converged"
                                        %folder)                    
                    
                    df0.loc[len(df0)] = [vac,
                                         cell,
                                         natoms,
                                         1/natoms,
                                         vr.final_energy,
                                         vr_ref.final_energy]
                    
                    
        df0.sort_values(['vacuum','N'],inplace=True)
        df0.to_excel(writer,'charge_0')
    

    ## modify dataframe for charged defects
    for q in [qi for qi in qs if qi != 'charge_0']:
        df = df0.copy(deep=True)
        
        for cell in myutils.listdironly(myutils.joinpath(args.path,q,'')):
            for vac in myutils.listdironly(myutils.joinpath(args.path,q,cell,'')):
                myLogger.info("parsing %s %s %s"%(q,cell,vac))

                folder = myutils.joinpath(args.path,q,cell,vac,'')

                if args.soc:  
                    folder = myutils.joinpath(folder,'dos','')
                    myLogger.info("parsing dos subdirectory")
                if os.path.exists(folder) and 'restart' in myutils.listdironly(folder):
                    folder = myutils.joinpath(folder,'restart','')
                    myLogger.info("parsing restart subdirectory")
                    
                vr_file = myutils.joinpath(folder,'vasprun.xml')
                
                if not os.path.exists(vr_file):
                    myLogger.warning("%s file does not exist"%vr_file)
                    
                else:
                    vr = Vasprun(vr_file)                    
                    if not vr.converged:
                        myLogger.warning("Vasp calculation in %s may not be converged"
                                        %folder)                       
                    df.loc[(df['vacuum'] == vac) & 
                           (df['supercell'] == cell),'E_def'] = vr.final_energy
        
        df.to_excel(writer, q)

    writer.save()
    
    myLogger.debug("Total time taken (s): %.2f"%(time.time()-time0))
             

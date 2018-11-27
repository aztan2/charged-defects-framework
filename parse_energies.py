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


if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description='Parse total energies from OUTCARs.')
    parser.add_argument('path',help='path to the directory containing all the output files')
    parser.add_argument('xlfile',help='excel filename to save the dataframe to')
    parser.add_argument('--logfile',help='logfile to save output to')
       
    ## read in the above arguments from command line
    args = parser.parse_args()
    
    ## set up logging
    if args.logfile:
        logging.basicConfig(filename=args.logfile,filemode='w',format='%(levelname)s:%(message)s', level=logging.DEBUG)    
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(logging.Formatter('%(levelname)s:%(message)s'))
        logging.getLogger('').addHandler(console)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

    time0 = time.time()
    
    writer = pd.ExcelWriter(args.path + args.xlfile)

    qs = listdironly(args.path)
    
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
        for cell in listdironly(args.path+'charge_0/'):
            for vac in listdironly(args.path+'charge_0/'+cell+'/'):   
                logging.info("parsing neutral %s %s"%(cell,vac))
                dir0 = args.path+'charge_0/'+cell+'/'+vac+'/'
                if not os.path.exists(dir0 + "vasprun.xml"):
                    logging.warning("vasprun.xml file does not exist in %s"%dir0)
                elif not os.path.exists(dir0 + "bulkref/vasprun.xml"):
                    logging.warning("vasprun.xml file does not exist in %s"%(dir0+"bulkref/"))
                else:
                    vr = Vasprun(dir0 + "vasprun.xml")
                    if not vr.converged:
                        logging.warning("Vasp calculation in %s may not be converged"%dir0)
                    vr_ref = Vasprun(dir0 + "bulkref/vasprun.xml")
                    if not vr.converged:
                        logging.warning("Vasp calculation in %s may not be converged"%(dir0+"bulkref/"))
                    natoms = np.sum(Poscar.from_file(dir0 + "bulkref/POSCAR").natoms)
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
#        qval = int(q.split('_')[-1])
        df = df0.copy(deep=True)
        for cell in listdironly(args.path+q+'/'):
            for vac in listdironly(args.path+q+'/'+cell+'/'):
                logging.info("parsing %s %s %s"%(q,cell,vac))
                dir1 = args.path+q+'/'+cell+'/'+vac+'/'
                if not os.path.exists(dir1 + "vasprun.xml"):
                    logging.warning("vasprun.xml file does not exist in %s"%dir1)
                else:
                    vr = Vasprun(dir1 + "vasprun.xml")
                    if not vr.converged:
                        logging.warning("Vasp calculation in %s may not be converged"%dir1)
                    df.loc[(df['vacuum'] == vac) & (df['supercell'] == cell),'E_def'] = vr.final_energy
        df.to_excel(writer, q)

    writer.save()
    
    logging.debug("Total time taken (s): %.2f"%(time.time()-time0))
               
#    df_1 = df_1.assign(new_col=df_1["E_def"] - df_1["E_bulk"])                
    
    
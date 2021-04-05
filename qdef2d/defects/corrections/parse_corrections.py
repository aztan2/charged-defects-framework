import os
import sys
import argparse
import numpy as np
import pandas as pd
from qdef2d import logging, osutils


def parse(dir_def,xlfile,soc=False,logfile=None):
    
    """
    Parse correction calculated by sxdefectalign2d and add to pandas dataframe.
    
    dir_def (str): path to the defect directory (also where the excel file is)
    xlfile (str): excel filename to read/save the dataframe from/to
    [optional] soc (bool): whether or not to look in soc(dos) subdirectory. Default=False.
    [optional] logfile (str): logfile to save output to
    
    """
    
    ## set up logging
    if logfile:
        myLogger = logging.setup_logging(logfile)
    else:
        myLogger = logging.setup_logging()
        

    ## load list of dataframes from sheets from excel file    
    df = pd.read_excel(os.path.join(dir_def,xlfile),sheet_name=None)
       
    for q in [qi for qi in df.keys() if qi != 'charge_0']:
        df[q]['E_corr'] = np.nan


    for q in osutils.listdironly(os.path.join(dir_def,'')):
        if q != 'charge_0':
            for cell in osutils.listdironly(os.path.join(dir_def,q,'')): 
                for vac in osutils.listdironly(os.path.join(dir_def,q,cell,'')): 
                    myLogger.info("parsing %s %s %s"%(q,cell,vac))
                    
                    folder = os.path.join(dir_def,q,cell,vac,'')
                    
                    if soc: 
                        myLogger.info("parsing dos subdirectory")
                        folder = os.path.join(folder,'dos','')     
                    if os.path.exists(folder) and 'restart' in osutils.listdironly(folder):
                        folder = os.path.join(folder,'restart','')
                        myLogger.info("parsing restart subdirectory")
    
                    if not os.path.exists(os.path.join(folder,'correction','correction')):
                        myLogger.warning("correction file does not exist")
                    else:
                        with open(os.path.join(folder,'correction','correction')) as f:
                            lines = f.readlines()
                            for line in lines:
                                if line[:21] == 'iso - periodic energy':
                                    E_corr = float(line.split()[-2])
                        df[q].loc[(df[q]['vacuum'] == vac) & 
                                  (df[q]['supercell'] == cell),'E_corr'] = E_corr


    ## write the updated excel file
    writer = pd.ExcelWriter(os.path.join(dir_def,xlfile))
    for q in df.keys():  
        df[q].to_excel(writer, q, index=False)
    writer.save()    
       
    
if __name__ == '__main__':
    
    
    ## this script can also be run directly from the command line
    parser = argparse.ArgumentParser(description='Parse correction calculated by sxdefectalign2d \
                                     and add to pandas dataframe.')
    parser.add_argument('dir_def',help='path to the defect directory (also where the excel file is)')
    parser.add_argument('xlfile',help='excel filename to read/save the dataframe from/to')
    parser.add_argument('--soc',help='whether or not to look in soc(dos) subdirectory',
                        default=False,action='store_true')
    parser.add_argument('--logfile',help='logfile to save output to')
       
    ## read in the above arguments from command line
    args = parser.parse_args()
    
    parse(args.dir_def, args.xlfile, args.soc, args.logfile) 


import os
import argparse
import logging
import numpy as np
import pandas as pd
import myutils


if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser(description='Parse correction calculated by sxdefectalign2d.')
    parser.add_argument('path',help='path to the directory containing the correction files')
    parser.add_argument('xlfile',help='excel filename to read/save the dataframe to')
    parser.add_argument('--soc',help='whether or not to look in soc(dos) subdirectory',default=False,action='store_true')
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


    ## load list of dataframes from sheets from excel file    
    df = pd.read_excel(myutils.joinpath(args.path,args.xlfile),sheet_name=None)
       
    for q in [qi for qi in df.keys() if qi != 'charge_0']:
        df[q]['E_corr'] = np.nan


    for q in myutils.listdironly(myutils.joinpath(args.path,'')):
        for cell in myutils.listdironly(myutils.joinpath(args.path,q,'')): 
            for vac in myutils.listdironly(myutils.joinpath(args.path,q,cell,'')): 
                logging.info("parsing %s %s %s"%(q,cell,vac))
                
                folder = myutils.joinpath(args.path,q,cell,vac,'')
                
                if args.soc: 
                    logging.info("parsing dos subdirectory")
                    folder = myutils.joinpath(folder,'dos','')     
                if os.path.exists(folder) and 'restart' in myutils.listdironly(folder):
                    folder = myutils.joinpath(folder,'restart','')
                    logging.info("parsing restart subdirectory")

                if not os.path.exists(myutils.joinpath(folder,'correction','correction')):
                    logging.warning("correction file does not exist")
                else:
                    with open(myutils.joinpath(folder,'correction','correction')) as f:
                        lines = f.readlines()
                        for line in lines:
                            if line[:21] == 'iso - periodic energy':
                                E_corr = float(line.split()[-2])
                    df[q].loc[(df[q]['vacuum'] == vac) & 
                              (df[q]['supercell'] == cell),'E_corr'] = E_corr
    

#    for cell in listdironly(joinpath(args.path,'')):
#        for vac in listdironly(joinpath(args.path,cell,'')): 
#            for q in listdironly(joinpath(args.path,cell,vac,'')): 
#                logging.info("parsing %s %s %s"%(cell,vac,q))
#                if not os.path.exists(joinpath(args.path,cell,vac,q,'correction')):
#                    logging.warning("correction file does not exist")
#                else:
#                    with open(joinpath(args.path,cell,vac,q,'correction')) as f:
#                        lines = f.readlines()
#                        for line in lines:
#                            if line[:21] == 'iso - periodic energy':
#                                E_corr = float(line.split()[-2])
#                    df[q].loc[(df[q]['vacuum'] == vac) & 
#                              (df[q]['supercell'] == cell),'E_corr'] = E_corr


    ## write the updated excel file
    writer = pd.ExcelWriter(myutils.joinpath(args.path,args.xlfile))
    for q in df.keys():  
        df[q].to_excel(writer, q)
    writer.save()    
       

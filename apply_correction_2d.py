import os
import sys
import json
import argparse
import myutils, gen_SPHInX_input_file, get_alignment_correction_2d


def main(args):
    
    ## define a main function callable from another python script if necessary
    
    parser = argparse.ArgumentParser(description='Apply sxdefectalign2d correction \
                                     to all charged defect calulations.')
    parser.add_argument('path_def', help='path to main defect directory')
    parser.add_argument('path_ref', help='path to the undefected reference directory')
    parser.add_argument('--eps_slab', type=float,
                        help='average slab dielectric constant (supply this or --file_dbentry)')
    parser.add_argument('--d_slab', type=float,
                        help='slab thickness (supply this or --file_dbentry)')
    parser.add_argument('--file_dbentry', help='path to the relevant dbentry.json file \
                        (supply this or --eps_slab and --d_slab)')
    parser.add_argument('--functional',  default='PBE',
                        help='functional that was used for this set of calculations')
    parser.add_argument('--encut', type=int, default=520, help='cutoff energy (eV)')
    parser.add_argument('--soc', default=False,action='store_true',
                        help='whether or not to look in soc(dos) subdirectory')
    parser.add_argument('--logfile', help='logfile to save output to')
       
    ## read in the above arguments from command line
    args = parser.parse_args()
    
    
    ## set up logging
    if args.logfile:
        myLogger = myutils.setup_logging(args.logfile)
    else:
        myLogger = myutils.setup_logging()

    
    ## if eps_slab and d_slab are directly provided
    if args.eps_slab and args.d_slab:
        eps_slab = args.eps_slab
        d_slab = args.d_slab
        myLogger.info("eps_slab: %.2f ; d_slab: %.2f"%(eps_slab,d_slab))

    elif args.file_dbentry:
        ## extract the eps_slab and d_slab from the relevant dbentry file
        if os.path.exists(args.file_dbentry):
            myLogger.info("Using slab properties from " + args.file_dbentry)
            with open(args.file_dbentry, 'r') as file:
                mater = json.loads(file.read())
                eps_slab = mater[args.functional]["unitcell"]["eps_ave"]
                d_slab = mater[args.functional]["unitcell"]["d_slab"]
                myLogger.info("eps_slab: %.2f ; d_slab: %.2f"%(eps_slab,d_slab))
        else:
            ## if can't find a dbentry file
            myLogger.info("Cannot find the file " + args.file_dbentry)
            
    else:
        myLogger.info("Insufficient information provided about the slab dielectric profile")


    qs = myutils.listdironly(args.path_def)
    for q in [qi for qi in qs if qi != 'charge_0']:      
        for cell in myutils.listdironly(myutils.joinpath(args.path_def,q,'')):
            for vac in myutils.listdironly(myutils.joinpath(args.path_def,q,cell,'')):
                folder = myutils.joinpath(args.path_def,q,cell,vac,'')
                folder_ref = myutils.joinpath(args.path_ref,'charge_0',cell,vac,'')

                if args.soc:  
                    folder = myutils.joinpath(folder,'dos','')
                    folder_ref = myutils.joinpath(folder_ref,'dos','')
                if os.path.exists(folder) and 'restart' in myutils.listdironly(folder):
                    folder = myutils.joinpath(folder,'restart','')
                
                ## check if defectproperty.json file is present in current directory
                if not os.path.exists(myutils.joinpath(folder,'defectproperty.json')):
                    myLogger.warning("defectproperty.json file does not exist; I can't proceed!")
                    
                else:
                    os.chdir(folder)
                    ## generate SPHInX input file
                    ## usage: gen_SPHInX_input_file.py [-h] eps slab_d
                    gen_SPHInX_input_file.main([str(eps_slab), str(d_slab)])
                    
                    if (os.path.exists(myutils.joinpath(folder,'LOCPOT')) and 
                        os.path.exists(myutils.joinpath(folder_ref,'LOCPOT'))):
                        
                        os.chdir(myutils.joinpath(folder,'correction'))
                        myLogger.info("applying correction to calculations in %s"%folder) 
                        
                        ## apply correction
                        ## usage: get_alignment_correction_2d.py [-h] vref vdef encut q
                        ##        [--threshold_slope THRESHOLD_SLOPE] [--threshold_C THRESHOLD_C]
                        ##        [--vfile VFILE] [--noplot] [--logfile LOGFILE]
                        get_alignment_correction_2d.main([myutils.joinpath(folder_ref,'LOCPOT'),
                                                          myutils.joinpath(folder,'LOCPOT'),
                                                          str(args.encut), q.split('_')[-1],
#                                                          '--allplots', 
                                                          '--logfile', 'getalign.log'])
    
                    elif not os.path.exists(myutils.joinpath(folder,'LOCPOT')):
                        myLogger.warning("Missing LOCPOT file in %s"%folder)
                        
                    elif not os.path.exists(myutils.joinpath(folder_ref,'LOCPOT')):
                        myLogger.warning("Missing LOCPOT file in %s"%folder_ref)    
                    

if __name__ == '__main__':

    
    main(sys.argv[1:])
    
    
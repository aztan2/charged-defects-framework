import os
import shutil
import sys
import argparse
from qdef2d import osutils, logging
from qdef2d.io.vasp import incar, kpoints
from qdef.defects import gen_defect_supercell


def main(dir_def_main,qs,cells,vacs,functional='PBE',kppa=400,write_bulkref=False):

    """ 
    Generate input files for defect calulations.
    
    Parameters
    ----------
    dir_def_main (str): path to the main defect directory
    qs (list of ints): list of charge states
    cells (list of strs): list of n1xn2xn3 supercell sizes
    vacs (list of ints): list of vacuum spacings
    [optional] functional (str): type of function: PBE(default)/SCAN+rVV10
    [optional] kppa (int): kpoint density per reciprocal atom. Default=400 pra.
    [optional] write_bulkref (str): write files for reference calculations? Default=False.
    
    """

    ## set up logging
    myLogger = logging.setup_logging()
    
    
    ## check if initdef file is present in dir_def_main ?
    if not write_bulkref:
        if osutils.check_file_exists(dir_def_main,"initdef") == True:
            for file in os.listdir(dir_def_main): 
                if file.startswith("initdef"):
                    file_initdef = file
        else:
            myLogger.warning("Stopping because of missing initdef file!")
            return
            
        
    ## check if POTCAR is present in dir_def_main ?
    if not osutils.check_file_exists(dir_def_main,"POTCAR"):
        myLogger.warning("Stopping because of missing POTCAR file!")
        return

    
    ## check if relevant POSCARs are present in dir_def_main ?
    for vac in vacs:
        if not osutils.check_file_exists(dir_def_main,"POSCAR_vac_%d"%vac) :
            myLogger.warning("Stopping because of missing POSCAR file(s)!")
            return            
    
    
    for q in qs:
        dir_q = os.path.join(dir_def_main,"charge_%d"%q)
        ## create and enter charge subdirectory
        if not os.path.exists(dir_q):
            os.makedirs(dir_q)
        os.chdir(dir_q)
        
        for cell in cells:
            dir_cell = os.path.join(dir_q,cell)
            ## create and enter supercell subdirectory
            if not os.path.exists(dir_cell):
                os.makedirs(dir_cell)
            os.chdir(dir_cell)
            
            for vac in vacs:
                dir_vac = os.path.join(dir_cell,"vac_%d"%vac)
                ## create and enter vacuum subdirectory
                if not os.path.exists(dir_vac):
                    os.makedirs(dir_vac)
                os.chdir(dir_vac)
                
                myLogger.info("generating input files in ",os.getcwd())
                
                ## cp POTCAR from dir_main
                shutil.copyfile(os.path.join(dir_def_main,"POTCAR"),
                                os.path.join(os.getcwd(),"POTCAR"))
                
                ## generate defect POSCAR
                if write_bulkref:
                    gen_defect_supercell([args.dir_main,'dummy.json',
                                                str(q),cell,str(vac),
                                                '--write_bulkref'])    ## TO UPDATE!!!
                else:
                    gen_defect_supercell([args.dir_main,
                                                (os.path.join(args.dir_main,file_initdef)),
                                                str(q),cell,str(vac)])     ## TO UPDATE!!!
                
                ## generate INCAR
                incar.gen_incar(q=q,functional=functional)
                
                ## generate KPOINTS
                kpoints.gen_kpts_uniform(kppa=kppa)
                
                ## generate submission script 
                ## with the default settings nodes=1, mem=2048, time=24:00:00
                if write_bulkref:
                    submit.gen_submit(jobname='ref_%s_%d'%(cell,vac))
                else:
                    submit.gen_submit()
                
    os.chdir(args.dir_main)


if __name__ == '__main__':

    
    ## this script can also be run directly from the command line    
    parser = argparse.ArgumentParser(description='Generate input files for defect calulations.')
    parser.add_argument('dir_main', help='absolute path to main defect directory')
    parser.add_argument('--qs', nargs='+', help='(required) charge states; \
                        list each charge state separated by a space', type=int)
    parser.add_argument('--cells', nargs='+', help='(required) supercell sizes; \
                        list each supercell size as n1xn2xn3 separated by a space')
    parser.add_argument('--vacs', nargs='+', help='(required) vacuum spacings; \
                        list each vacuum spacing separated by a space', type=int)
    parser.add_argument('--functional', default='PBE',
                        help='type of function: PBE(default)/SCAN+rVV10')
#    parser.add_argument('--soc',default=False,action='store_true',
#                        help='turn on spin-orbit coupling')
    parser.add_argument('--kppa', type=int, help='kpt density (pra)', default=440)
    parser.add_argument('--write_bulkref',help='to write bulkref?',
                        default=False,action='store_true')
      
    ## parse the given arguments
    args = parser.parse_args(args)
    
    
    main(sys.argv[1:])
    

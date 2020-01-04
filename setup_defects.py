import os
import shutil
import sys
import argparse
import gen_defect_supercells, gen_incar, gen_kpts_grid, gen_submit


def check_file_exists(directory,filename):
    
    ## check if a file is present in a given directory
    
    files = [f for f in os.listdir(directory) 
                if f.startswith(filename)]
    if len(files) < 1:
        raise RuntimeError ("can't find %s file!"%filename)
        sys.exit(1)
    elif len(files) > 1:
        raise RuntimeError ("more than 1 %s file found"%filename)
        sys.exit(1)
    else:
        return True


def main(args):
    
    ## define a main function callable from another python script if necessary
    
    parser = argparse.ArgumentParser(description='Generate all input files \
                                     for defect calulations.')
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
 
    
    ## check if initdef file is present in dir_main ?
    if not args.write_bulkref:
        if check_file_exists(args.dir_main,"initdef") == True:
            for file in os.listdir(args.dir_main): 
                if file.startswith("initdef"):
                    file_initdef = file
        
    ## check if POTCAR is present in dir_main ?
    check_file_exists(args.dir_main,"POTCAR")  

    ## check if relevant POSCARs are present in dir_main ?
    for vac in args.vacs:
        check_file_exists(args.dir_main,"POSCAR_vac_%d"%vac)        
    
    
    for q in args.qs:
        dir_q = os.path.join(args.dir_main,"charge_%d"%q)
        ## create and enter charge subdirectory
        if not os.path.exists(dir_q):
            os.makedirs(dir_q)
        os.chdir(dir_q)
        
        for cell in args.cells:
            dir_cell = os.path.join(dir_q,cell)
            ## create and enter supercell subdirectory
            if not os.path.exists(dir_cell):
                os.makedirs(dir_cell)
            os.chdir(dir_cell)
            
            for vac in args.vacs:
                dir_vac = os.path.join(dir_cell,"vac_%d"%vac)
                ## create and enter vacuum subdirectory
                if not os.path.exists(dir_vac):
                    os.makedirs(dir_vac)
                os.chdir(dir_vac)
                print (os.getcwd())
                
                ## cp POTCAR from dir_main
                shutil.copyfile(os.path.join(args.dir_main,"POTCAR"),
                                os.path.join(os.getcwd(),"POTCAR"))
                
                ## generate defect POSCAR
                if args.write_bulkref:
                    gen_defect_supercells.main([args.dir_main,'dummy.json',
                                                str(q),cell,str(vac),
                                                '--write_bulkref'])
                else:
                    gen_defect_supercells.main([args.dir_main,
                                                (os.path.join(args.dir_main,file_initdef)),
                                                str(q),cell,str(vac)])
                
                ## generate INCAR
                gen_incar.main(['--q',str(q),'--functional',args.functional])
                
                ## generate KPOINTS
                gen_kpts_grid.main(['--kppa',str(args.kppa)])
                
                ## generate submission script 
                ## with the default settings nodes=1, mem=2048, time=2-00:00:00
                if args.write_bulkref:
                    gen_submit.main(['--jobname',"ref_%s_%d"%(cell,vac)])
                else:
                    gen_submit.main([])
                
    os.chdir(args.dir_main)


if __name__ == '__main__':


#    cmd = 'C:/Users/Anne/Desktop/UF/research/charged_defects/testing/bulkref/GGA/mag/ \
#            --q 0 --cell 3x3x1 4x4x1 --vac 15 20 --write_bulkref'
    
    cmd = 'C:/Users/Anne/Desktop/UF/research/charged_defects/testing/Sevac_/GGA/mag/ \
            --q 0 -1 1 --cell 3x3x1 4x4x1 --vac 15 20'
            
    main(cmd.split())
    
#    main(sys.argv[1:])

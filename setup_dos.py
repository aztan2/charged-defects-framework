import os
import shutil
import sys
import argparse
import gen_incar


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
    parser.add_argument('--soc',default=False,action='store_true',
                        help='turn on spin-orbit coupling')

      
    ## parse the given arguments
    args = parser.parse_args(args)        
    
    
    for q in args.qs:
        dir_q = os.path.join(args.dir_main,"charge_%d"%q)
        ## enter charge subdirectory
        os.chdir(dir_q)
        
        for cell in args.cells:
            dir_cell = os.path.join(dir_q,cell)
            ## enter supercell subdirectory
            os.chdir(dir_cell)
            
            for vac in args.vacs:
                dir_vac = os.path.join(dir_cell,"vac_%d"%vac)
                ## enter vacuum subdirectory
                os.chdir(dir_vac)
                ## create and enter dos subdirectory
                dir_dos = os.path.join(dir_vac,"dos")
                if not os.path.exists(dir_dos):
                    os.makedirs(dir_dos)
                os.chdir(dir_dos)
                print (os.getcwd())
                
                ## cp CONTCAR, KPOINTS, POTCAR, submitVASP.sh, defectproperty.json
                ## from converged relaxation run
                shutil.copyfile(os.path.join(dir_vac,"CONTCAR"),
                                os.path.join(os.getcwd(),"POSCAR"))
                for file in ["KPOINTS","POTCAR","submitVASP.sh","defectproperty.json"]:
                    shutil.copyfile(os.path.join(dir_vac,file),
                                    os.path.join(os.getcwd(),file))
                
                ## generate INCAR for dos calculation
                if args.soc:
                    gen_incar.main(['--q',str(q),'--runtype','dos','--functional',args.functional,'--soc'])
                else:
                    gen_incar.main(['--q',str(q),'--runtype','dos','--functional',args.functional])
                
                
    os.chdir(args.dir_main)


if __name__ == '__main__':

    
    main(sys.argv[1:])
    
    
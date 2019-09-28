import os
import json
import argparse
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pymatgen
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.io.vasp.inputs import Poscar



def sbatch_cmds(s,jobname,nodes,mem,time,qos):      
    
    s += '#SBATCH --job-name=%s\n'%jobname
    s += '#SBATCH -o out_%j\n'
    s += '#SBATCH -e err_%j\n' 
    s += '#SBATCH --mail-type=END,FAIL  # Mail events (NONE, BEGIN, END, FAIL, ALL)\n'
    s += '#SBATCH --mail-user=annemarietan@ufl.edu\n'
    s += '#SBATCH --partition=hpg2-compute\n'
    s += '#SBATCH --qos=%s\n'%qos
    s += '#SBATCH --ntasks=%d\n'%(nodes*32)
    s += '#SBATCH --ntasks-per-socket=16\n'
    s += '#SBATCH --ntasks-per-node=32\n'
    s += '#SBATCH --nodes=%d\n'%nodes
    s += '#SBATCH --cpus-per-task=1\n'
    s += '#SBATCH --mem-per-cpu=%dmb\n'%(mem)
    s += '#SBATCH -t %s\n'%time
    s += '#SBATCH --distribution=cyclic:cyclic\n\n'
    
    return s


def load_modules(s,vasp="noz"):      

    s += 'module purge\n'
    if vasp == "noz":
        s += 'module load intel/2018\n'
        s += 'module load openmpi/3.1.0\n'
        s += 'srun --mpi=pmix_v1 /home/annemarietan/vasp_noz > job.log\n\n'
    elif vasp == "ncl_noz":
        s += 'module load intel/2018\n'
        s += 'module load openmpi/3.1.0\n'
        s += 'srun --mpi=pmix_v1 /home/annemarietan/vasp_ncl_noz > job.log\n\n'
    elif vasp == "std":
        s += 'module load intel/2018\n'
        s += 'module load openmpi/3.0.0\n'
        s += 'module load vasp/5.4.4\n'
        s += 'srun --mpi=pmix_v1 vasp_std > job.log\n\n'
    else:
        raise ValueError ("unrecognized vasp version!")
    
    return s



if __name__ == '__main__':
 

    parser = argparse.ArgumentParser(description='Generate submit script')
    parser.add_argument('--jobname',help='jobname',default='jobby')
    parser.add_argument('--queue',help='queue',default='hennig')
    parser.add_argument('--nodes',type=int,help='number of nodes',default=1)
    parser.add_argument('--mem',type=int,help='memory per node',default=2048)
    parser.add_argument('--time',help='time requested (d-hh:mm:ss)',default='2-00:00:00')
    parser.add_argument('--vasp',help='vasp executable (noz/ncl_noz/std)',default='noz')
      
    ## read in the above arguments from command line
    args = parser.parse_args()
    
    
    ## the bash script already put us in the appropriate subdirectory
    dir_sub = os.getcwd()
 
    if os.path.exists(os.path.join(dir_sub,"defectproperty.json")):
        with open(os.path.join(dir_sub,"defectproperty.json"), 'r') as file:
            defprop = json.loads(file.read())
        jobname = ("_".join(defprop["defect_type"])
                    +"_"+str(defprop["charge"])
                    +"_"+"x".join([str(m) for m in defprop["supercell"]])
                    +"_"+str(defprop["vacuum"]))
    else:
        jobname = args.jobname
    print (jobname)
    
    s = '#!/bin/bash\n'
    s = sbatch_cmds(s,jobname,args.nodes,args.mem,args.time,args.queue)
    s += 'cd $SLURM_SUBMIT_DIR\n\n'
    s = load_modules(s,args.vasp)
    s += 'echo \'Done.\''
    
    with open(os.path.join(dir_sub,"submitVASP.sh"),'w') as f:
        f.write(s)

    
#    dir_main = "mp-2815_MoS2/test/GGA/mag/"
#    dir_main = "Y:/MoS2/monolayer_Nbsub/SCAN_vdW/mag/"
#    dir_main = "Y:/WSe2/monolayer/SCAN_vdW/mag/"
       
#    for q in [0,-1]:
#        for cell in [(4,4)]:
#            for vac in [15,20]:
#                
#                dir_sub = dir_main+"charge_%d/%dx%dx1/vac_%d/"%(q,cell[0],cell[1],vac)
#               
#                with open(dir_sub+"defectproperty.json", 'r') as file:
#                    defprop = json.loads(file.read())
#                jobname = ("_".join(defprop["defect_type"])
#                            +"_"+str(defprop["charge"])
#                            +"_"+"x".join([str(m) for m in defprop["supercell"]])
#                            +"_"+str(defprop["vacuum"]))
#                print (jobname)
#                
#                s = '#!/bin/bash\n'
#                s = sbatch_cmds(s,jobname=jobname,nodes=2)
#                s += 'cd $SLURM_SUBMIT_DIR\n\n'
#                s = load_modules(s)
#                s += 'echo \'Done.\''
#                
#                with open(dir_sub+"submitVASP.sh",'w') as f:
#                    f.write(s)
                    
                    
#    for vac in [10,12,14,15,16,18,20]:
#        
#        dir_sub = dir_main+"vac_%d/"%(vac)
#       
#        jobname = "WSe2_vac"+str(vac)
#        print (jobname)
#        
#        s = '#!/bin/bash\n'
#        s = sbatch_cmds(s,jobname=jobname,nodes=1)
#        s += 'cd $SLURM_SUBMIT_DIR\n\n'
#        s = load_modules(s)
#        s += 'echo \'Done.\''
#        
#        with open(dir_sub+"submitVASP.sh",'w') as f:
#            f.write(s)
            
            
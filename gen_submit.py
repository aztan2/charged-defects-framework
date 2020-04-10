import sys
import os
import json
import argparse


def sbatch_cmds(s,jobname,nodes,mem,time,qos,email=None):
    
    s += '#SBATCH --job-name=%s\n'%jobname
    s += '#SBATCH -o out_%j\n'
    s += '#SBATCH -e err_%j\n'
    if email is not None:
        s += '#SBATCH --mail-type=END,FAIL  # Mail events (NONE, BEGIN, END, FAIL, ALL)\n'
        s += '#SBATCH --mail-user=%s\n'%email
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
    if vasp == "noz_intel2019":
        s += 'module load intel/2019.1.144\n'
        s += 'module load openmpi/4.0.1\n'
        s += 'srun --mpi=pmix_v3 /home/annemarietan/vasp5.4.4_intel2019/vasp_noz > job.log\n\n'
    elif vasp == "ncl_intel2019":
        s += 'module load intel/2019.1.144\n'
        s += 'module load openmpi/4.0.1\n'
        s += 'srun --mpi=pmix_v3 /home/annemarietan/vasp5.4.4_intel2019/vasp_ncl > job.log\n\n'
    elif vasp == "noz":
        s += 'module load intel/2018\n'
        s += 'module load openmpi/3.1.0\n'
        s += 'srun --mpi=pmix_v1 /home/annemarietan/vasp5.4.4_intel2018/vasp_noz > job.log\n\n'
    elif vasp == "ncl_noz":
        s += 'module load intel/2018\n'
        s += 'module load openmpi/3.1.0\n'
        s += 'srun --mpi=pmix_v1 /home/annemarietan/vasp5.4.4_intel2018/vasp_ncl_noz > job.log\n\n'
    elif vasp == "std":
        s += 'module load intel/2019.1.144\n'
        s += 'module load openmpi/4.0.0\n'
        s += 'module load vasp/5.4.4\n'
        s += 'srun --mpi=pmix_v3 vasp_std > job.log\n\n'
#        s += 'module load intel/2018\n'
#        s += 'module load openmpi/3.0.0\n'
#        s += 'module load vasp/5.4.4\n'
#        s += 'srun --mpi=pmix_v1 vasp_std > job.log\n\n'
    else:
        raise ValueError ("unrecognized vasp version!")
    
    return s


def main(args):
    
    ## define a main function callable from another python script
    
    parser = argparse.ArgumentParser(description='Generate submit script')
    parser.add_argument('--jobname',help='jobname',default='jobby')
    parser.add_argument('--queue',help='queue',default='hennig')
    parser.add_argument('--nodes',type=int,help='number of nodes',default=1)
    parser.add_argument('--mem',type=int,help='memory per node',default=2048)
    parser.add_argument('--time',help='time requested (d-hh:mm:ss)',default='2-00:00:00')
    parser.add_argument('--vasp',help='vasp executable (noz/ncl_noz/std)',default='noz_intel2019')
    parser.add_argument('--email',help='user email',default='annemarietan@ufl.edu')
      
    ## read in the above arguments from command line
    args = parser.parse_args(args)
    
    
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
    s = sbatch_cmds(s,jobname,args.nodes,args.mem,args.time,args.queue,args.email)
    s += 'cd $SLURM_SUBMIT_DIR\n\n'
    s = load_modules(s,args.vasp)
    s += 'echo \'Done.\''
    
    with open(os.path.join(dir_sub,"submitVASP.sh"),'w') as f:
        f.write(s)


if __name__ == '__main__':
 
    main(sys.argv[1:])  
            
            

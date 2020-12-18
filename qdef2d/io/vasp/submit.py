import os
import json


def sbatch_cmds(s,jobname,qos,nodes,cpus,mem,time,email=None):
    
    s += '#SBATCH --job-name=%s\n'%jobname
    s += '#SBATCH -o out_%j\n'
    s += '#SBATCH -e err_%j\n'
    if email is not None:
        s += '#SBATCH --mail-type=END,FAIL  # Mail events (NONE, BEGIN, END, FAIL, ALL)\n'
        s += '#SBATCH --mail-user=%s\n'%email
    s += '#SBATCH --partition=hpg2-compute\n'
    s += '#SBATCH --qos=%s\n'%qos
    s += '#SBATCH --ntasks=%d\n'%cpus
#    s += '#SBATCH --ntasks-per-socket=16\n'
#    s += '#SBATCH --ntasks-per-node=32\n'
    s += '#SBATCH --nodes=%d\n'%nodes
    s += '#SBATCH --cpus-per-task=1\n'
    s += '#SBATCH --mem-per-cpu=%dmb\n'%(mem)
    s += '#SBATCH -t %s\n'%time
    s += '#SBATCH --distribution=cyclic:cyclic\n\n'
    
    return s


def load_modules(s,vasp="noz_intel2019"):      

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


def gen_submit(jobname='jobby',queue='hennig',nodes=1,cpus=None,mem=2048,time='24:00:00',
               vasp='noz_intel2019',email='annemarietan@ufl.edu'):
    
    """ 
    Generate submission script for hipergator.
    
    Parameters
    ----------
    [optional] jobname (str): job name. Default='jobby'.
    [optional] queue (str): queue. Default='hennig'.
    [optional] nodes (int): number of nodes requested. Default=1.
    [optional] cpus (int): number of cpus requested. Default=nodes*32.
    [optional] mem (int): requested memory per node. Default=2048(MB). 
    [optional] time (str): time requested (d-hh:mm:ss). Default='24:00:00'.
    [optional] vasp (str): type of vasp executable to use, e.g. noz/ncl_noz/std. Default='noz_intel2019'.
    [optional] email (str): user email. Default='annemarietan@ufl.edu'. PLEASE CHANGE ACCORDINGLY!!
       
    """
    
    dir_sub = os.getcwd()
        
    if os.path.exists(os.path.join(dir_sub,"defectproperty.json")):
        with open(os.path.join(dir_sub,"defectproperty.json"), 'r') as file:
            defprop = json.loads(file.read())
        jobname = ("_".join(defprop["defect_type"])
                    +"_"+str(defprop["charge"])
                    +"_"+"x".join([str(m) for m in defprop["supercell"]])
                    +"_"+str(defprop["vacuum"]))
    else:
        jobname = jobname
    print (jobname)
    
    if cpus == None: cpus = nodes*32
    
    s = '#!/bin/bash\n'
    s = sbatch_cmds(s,jobname,queue,int(nodes),int(cpus),int(mem),time,email)
    s += 'cd $SLURM_SUBMIT_DIR\n\n'
    s = load_modules(s,vasp)
    s += 'echo \'Done.\''
    
    with open(os.path.join(dir_sub,"submitVASP.sh"),'w') as f:
        f.write(s) 

        
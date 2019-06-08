import os
import argparse
import numpy as np
from pymatgen.io.vasp.inputs import Poscar,Potcar,Incar


class IncarSettings:
    
    def __init__(self,func,runtype,charge,poscar,potcar,wavecar,chgcar):
        
        self.func = func
        self.runtype = runtype
        self.charge = charge
        self.poscar = poscar
        self.potcar = potcar
        self.chgcar = chgcar
        self.wavecar = wavecar
        self.params = {}
        
        
    def setparams(self):
        
        self.functional()
        self.startup()
        self.elrelax()
        self.ionicrelax()
        self.smearing()
        self.mag()
        self.parallel()
        self.output()
        
        return
        

    def functional(self):
    
        if "SCAN" in self.func:
            self.params.update({"METAGGA": "SCAN",
                                "LASPH": True,
                                "ADDGRID": True
                           })
        if "rVV10" in self.func:
            self.params.update({"BPARAM": 15.7,
                                "LUSE_VDW": True
                                })
        if "optB88-vdW" in self.func:
            self.params.update({"GGA": "BO",
                                "PARAM1": 0.1833333333,
                                "PARAM2": 0.2200000000,
                                "LUSE_VDW": True,
                                "AGGAC": 0.0000
                                })
        
        return
    
    
    def startup(self,prec='Accurate',algo='Fast',lreal='Auto',
                istart=None,icharg=None,isym=0,nelect=None,encut=520):

        if "SCAN" in self.func:
            algo = 'All'
            
        if self.runtype in ("dos","bands"):
            if self.wavecar:
                icharg = 10  ## Read charge density from WAVECAR and keep fixed
            elif self.chgcar:
                icharg = 11  ## Read charge density from CHGCAR and keep fixed
                
        nelect = int(np.dot([pc.nelectrons for pc in self.potcar],
                            self.poscar.natoms) - self.charge)
        
        params = {"PREC": prec,
                  "ALGO": algo,  ## switch to 'All' if using SCAN
                  "LREAL": lreal,
                  "ISTART": istart,
                  "ICHARG": icharg,  ## set to 10 or 11 for dos or bands calcs 
                  "ISYM": isym,  ## set to 0 to turn off symmetry
                  "NELECT": nelect,  ## specify for charged systems
                  "ENCUT": encut
                  }
        
        self.params.update(params)
        
        return
    
        
    def elrelax(self,nelm=120,ediff=1e-6):
        
        params = {"NELM": nelm, ## max number of electronic steps per SC loop
                  "EDIFF": ediff  ## accuracy for electronic minimization
                  }
        
        self.params.update(params)
        
        return
    
    
    def ionicrelax(self,isif=2,ibrion=2,nsw=100):
        
        if self.runtype in ("dos","bands"):
            ibrion = -1
            nsw = 0

        params = {"ISIF": isif,  ## switch to 3 to relax cell vectors
                  "IBRION": ibrion,  ## CG relaxation; switch to -1 for static calc
                  "NSW": nsw  ## max number of ionic steps
                  }
        
        self.params.update(params)
        
        return
        
     
    def smearing(self,ismear=1,sigma=0.1,emin=None,emax=None,nedos=None):
        
        if self.runtype == "dos":
            ismear = 0  ## tetrahedron smearing (-5) sometimes gives problems?
            sigma = 0.01
            emin = -20
            emax = 10
            nedos = 3001
                     
        params = {"ISMEAR": ismear,  ## MP smearing; switch to -5 or 0 for dos calc
                  "SIGMA": sigma,  ## smearing width
                  "EMIN": emin,  ## specify for dos calc
                  "EMAX": emax,  ## specify for dos calc
                  "NEDOS": nedos  ## specify for dos calc
                  }
        
        self.params.update(params)
        
        return
    
    
    def mag(self,ispin=2):

        if ispin == 2:
            high_magmom, low_magmom = 0,0
            for pc,atomnum in zip(potcar,poscar.natoms):
                if 'd' in [orb[1] for orb in pc.electron_configuration]:
                    high_magmom += atomnum
                else:
                    low_magmom += atomnum
                    
            params = {"ISPIN": ispin,  ## default is spin-polarized
                      "MAGMOM": "%d*5.0 %d*0.6"%(high_magmom,low_magmom)
                      }
        
        self.params.update(params)
        
        return
        
        
    def parallel(self,lplane=True,npar=4,kpar=2):
        
        params = {"LPLANE": True,  ## parallelize over planewaves
                  "NPAR": npar,  ## number of bands treated in parallel (~sqrt #nodes)
                  "KPAR": kpar   ## number of k-points treated in parallel 
                  }
        
        self.params.update(params)
        
        return
    
    
    def output(self,lwave=True,lcharg=True,lmaxmix=4,lorbit=11,lvtot=True,lvhar=True):
        
        params = {"LWAVE": lwave,  ## write WAVECAR
                  "LCHARG": lcharg,   ## write CHGCAR
                  "LMAXMIX": lmaxmix,  ##  controls up to which l-quantum number are written to CHGCAR
                  "LORBIT": lorbit,  ## write PROCAR (projected wavefunctions)
                  "LVTOT": lvtot,   ## write LOCPOT
                  "LVHAR": lvhar
                  }
        
        self.params.update(params)
        
        return
        
        
    def stripNone(self):
        
        for param,val in list(self.params.items()):
            if val is None:
                del self.params[param]
      
        return


if __name__ == '__main__':
    
 
    parser = argparse.ArgumentParser(description='Generate INCAR')
    parser.add_argument('q',type=int,help='charge')
    parser.add_argument('--runtype',help='type of calculation: relax(default)/dos/bands',default='relax')
    parser.add_argument('--functional',help='type of function: PBE(default)/SCAN+rVV10',default='PBE')
      
    ## read in the above arguments from command line
    args = parser.parse_args()
    

    ## the bash script already put us in the appropriate subdirectory
    dir_sub = os.getcwd()
                
    if args.runtype == "dos":
        dir_sub = dir_sub+"dos/"
    
    poscar = Poscar.from_file(os.path.join(dir_sub,"POSCAR"))
    potcar = Potcar.from_file(os.path.join(dir_sub,"POTCAR"))
    wavecar = os.path.exists(os.path.join(dir_sub,"WAVECAR"))  ## just a boolean (does file exist)
    chgcar = os.path.exists(os.path.join(dir_sub,"CHGCAR"))  ## just a boolean (does file exist)

    inc = IncarSettings(args.functional,args.runtype,args.q,poscar,potcar,wavecar,chgcar)
    inc.setparams()
    inc.stripNone()       
    
    with open(os.path.join(dir_sub,"INCAR"),'w') as f:
        f.write(Incar.get_string(inc.params,sort_keys=False))
    
    
#    func = "PBE"
##    func = ["SCAN","rVV10"]
#    if type(func) is not list:
#        func = [func]
##    runtype = "relax"
##    runtype = "dos"
#    runtype = "bands"
#    
##    dir_main = "mp-2815_MoS2/"
##    dir_main = "Y:/MoS2/monolayer_Nbsub/SCAN_vdW/mag/"
#    dir_main = "Y:/WSe2/monolayer/GGA/mag/"
#    
#    for q in [-2]:
#        for cell in [(4,4),(5,5)]:
#            for vac in [15,20]:
#                
#                dir_sub = dir_main+"charge_%d/%dx%dx1/vac_%d/"%(q,cell[0],cell[1],vac)
#                if runtype == "dos":
#                    dir_sub = dir_sub+"dos/"
#                
#                poscar = Poscar.from_file(dir_sub+"POSCAR")
#                potcar = Potcar.from_file(dir_sub+"POTCAR")
#                wavecar = os.path.exists(dir_sub+"WAVECAR")  ## just a boolean (does file exist)
#                chgcar = os.path.exists(dir_sub+"CHGCAR")  ## just a boolean (does file exist)
#
#                inc = IncarSettings(func,runtype,q,poscar,potcar,wavecar,chgcar)
#                inc.setparams()
#                inc.stripNone()       
#                
#                with open(dir_sub+"INCAR",'w') as f:
#                    f.write(Incar.get_string(inc.params,sort_keys=False))


#    ## unitcell
#    for vac in [10,12,14,15,16,18,20]:
#        
#        q = 0
#        
#        dir_sub = dir_main+"vac_%d/"%(vac)
#        if runtype == "dos":
#            dir_sub = dir_sub+"dos/"
#        if runtype == "bands":
#            dir_sub = dir_sub+"bands/"
#        
#        poscar = Poscar.from_file(dir_sub+"POSCAR")
#        potcar = Potcar.from_file(dir_sub+"POTCAR")
#        wavecar = os.path.exists(dir_sub+"WAVECAR")  ## just a boolean (does file exist)
#        chgcar = os.path.exists(dir_sub+"CHGCAR")  ## just a boolean (does file exist)
#
#        inc = IncarSettings(func,runtype,q,poscar,potcar,wavecar,chgcar)
#        inc.setparams()
#        inc.ionicrelax(isif=3)
#        inc.parallel(lplane=True,npar=None,kpar=None)
#        inc.stripNone()       
#        
#        with open(dir_sub+"INCAR",'w') as f:
#            f.write(Incar.get_string(inc.params,sort_keys=False))

      
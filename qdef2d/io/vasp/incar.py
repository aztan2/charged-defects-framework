import os
import numpy as np
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pymatgen
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
            else:
                if self.runtype == "bands":
                    raise FileNotFoundError('No CHGCAR or WAVECAR provided '\
                                            'for a bandstructure calculation. '\
                                            'If you wish to proceed, you will need '\
                                            'to generate a special KPOINTS file '\
                                            'with a uniform mesh at the top.')
                
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
    
        
    def elrelax(self,nelm=120,ediff=1e-5):
        
        params = {"NELM": nelm, ## max number of electronic steps per SC loop
                  "EDIFF": ediff  ## accuracy for electronic minimization
                  }
        
        self.params.update(params)
        
        return
    
    
    def ionicrelax(self,isif=2,ibrion=2,nsw=100):
        
        if self.runtype in ("dos","bands"):
            ibrion = -1
            nsw = 0
            
        if self.runtype == "dielectric":
            ## get the ionic contribution to the macroscopic
            ## dielectric tensor using finite differences
            ibrion = 6  

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
    
    
    def mag(self,ispin=2,ncl=False):

        if ispin == 2:
            high_magmom, low_magmom = 0,0
            for pc,atomnum in zip(self.potcar,self.poscar.natoms):
                if 'd' in [orb[1] for orb in pc.electron_configuration]:
                    high_magmom += atomnum
                else:
                    low_magmom += atomnum
                    
            if ncl:  ## define xyz magmom for each atom    
                params = {"ISPIN": ispin,  ## default is spin-polarized
                          "MAGMOM": "%d*3.0 %d*0.35"%(3*high_magmom,3*low_magmom)
                          }
            else:        
                params = {"ISPIN": ispin,  ## default is spin-polarized
                          "MAGMOM": "%d*5.0 %d*0.6"%(high_magmom,low_magmom)
                          }
        
        self.params.update(params)
        
        return


    def soc(self,icharg=None,saxis=[1,1,1]):
        
        if self.chgcar:
            icharg = 1  ## Read charge density from CHGCAR (don't keep fixed)
            saxis = [0,0,1]

        params = {"ICHARG": icharg,
                  "LSORBIT": True,  ## turn on spin-orbit coupling
                  "LNONCOLLINEAR": True,  ## turn on noncollinear calculation
                  "SAXIS": "%d %d %d"%(saxis[0],saxis[1],saxis[2]),  ## specify reference axis
                  "GGA_COMPAT": False   ## do not restore full lattice symmetry
                  }
        
        if self.chgcar:
            params.update({"MAGMOM": None})  ## reads initial magmom from CHGCAR

        self.params.update(params)
        
        return
    
    
    def dielectric(self):
        
        if "SCAN" in self.func:
            ## get the macroscopic dielectric tensor using finite field method  
            params = {"LCALCEPS": True,
                      "LPEAD": True
                      }  
        else:
            ## get the macroscopic dielectric tensor using DFPT 
            params = {"LEPSILON": True,
                      "LPEAD": True
                      }  
            
        self.params.update(params)
        
        return
        
        
    def parallel(self,lplane=True,ncore=8):
        
        params = {"LPLANE": True,  ## parallelize over planewaves
                  "NCORE": ncore ## number of cores that work on an individual orbital
#                  "NPAR": npar,  ## number of bands treated in parallel (~sqrt #nodes)
#                  "KPAR": kpar   ## number of k-points treated in parallel 
                  }
        
        self.params.update(params)
        
        return
    
    
    def output(self,lwave=False,lcharg=False,lmaxmix=4,lorbit=11,lvtot=True,lvhar=True):
        
        params = {"LWAVE": lwave,  ## write WAVECAR
                  "LCHARG": lcharg,   ## write CHGCAR
                  "LMAXMIX": lmaxmix,  ##  controls up to which l-quantum number are written to CHGCAR
                  "LORBIT": lorbit,  ## write PROCAR (projected wavefunctions)
#                  "LVTOT": lvtot,   ## write LOCPOT
                  "LVHAR": lvhar   ## write LOCPOT
                  }
        
        self.params.update(params)
        
        return
        
        
    def stripNone(self):
        
        for param,val in list(self.params.items()):
            if val is None:
                del self.params[param]
      
        return
    
    
def generate(q=0,runtype='relax',functional='PBE',soc=False,relaxcell=False):

    """ 
    Generate INCAR file.
    
    Parameters
    ----------
    [optional] q (int): net charge of system. Default=0.
    [optional] runtype (str): type of calculation: relax(default)/dos/bands/dielectric
    [optional] functional (str): type of functional: PBE(default)/SCAN+rVV10
    [optional] soc (bool): to include spin-orbit coupling? Default=False.
    [optional] relaxcell (bool): to relax cell lattice parameters? Default=Fase.
    
    """
    
    dir_sub = os.getcwd()
    
    poscar = Poscar.from_file(os.path.join(dir_sub,"POSCAR"))
    potcar = Potcar.from_file(os.path.join(dir_sub,"POTCAR"))
    wavecar = os.path.exists(os.path.join(dir_sub,"WAVECAR"))  ## just a boolean (does file exist)
    chgcar = os.path.exists(os.path.join(dir_sub,"CHGCAR"))  ## just a boolean (does file exist)

    inc = IncarSettings(functional.split("+"),runtype,int(q),poscar,potcar,wavecar,chgcar)
    inc.setparams()
    if soc:
        inc.mag(ncl=True)
        inc.soc()
    if relaxcell:
        inc.ionicrelax(isif=3)
        inc.parallel(ncore=1)
        inc.output(lcharg=True)
    if runtype == 'dielectric':
        inc.dielectric()
        inc.startup(isym=None)
        inc.parallel(ncore=None)
        inc.output(lvtot=None,lvhar=None)
    inc.stripNone()      
    
    with open(os.path.join(dir_sub,"INCAR"),'w') as f:
        f.write(Incar.get_string(inc.params,sort_keys=False))
   

import re
import os
import sys
import json
import argparse
from pymatgen.io.vasp.outputs import Poscar, Vasprun
import myutils


class DatabaseEntry(object):
    
    """
    This class contains the functionalities related to
    creating, manipulating, and reading simple database entries
    
    """
    
    def __init__(self, system, monolayer, dir_dft, dir_db, dbentry, funcs, myLogger):
        
        self.system = system
        self.monolayer = monolayer
        self.dir_dft = dir_dft
        self.dir_db = dir_db
        self.dbentry = dbentry
        self.funcs = funcs
        self.log = myLogger
        
        return


    def create_entry_from_vasp(self, eps_slab=None, d_slab=None):
        
        """
        create a database entry from data obtained directly from vasp output files
        with the format 
        {
            functional: {
                mu: xxx, 
                vac_x: {
                    VBM: xxx
                }
            }
        } 
        
        """
        
        if self.dbentry:
            ## load existing json file
            mater = self.load_from_json()
        else:
            ## or initialize empty dict
            mater = {}
            self.dbentry = "%s.json"%self.system
        
        for i,func in enumerate(self.funcs):
            
            if func not in mater:
                mater[func] = {}
        
            dir_func = os.path.join(self.dir_dft,func.split('+')[0])
            if "mag" in myutils.listdironly(dir_func):
                dir_func =  os.path.join(dir_func,"mag")
                
            if self.monolayer:
                ## for a monolayer system, enter every vacuum subdirectory
                ## and extract total energy, Band gap, VBM
                for vac in myutils.listdironly(dir_func):
                    
                    if vac[0:3] == "vac":
                        dir_vac = os.path.join(dir_func,vac)
                        if func.split('+')[-1] == "soc":
                             dir_vac = os.path.join(dir_vac,"soc")

                        vrfile = os.path.join(dir_vac,'vasprun.xml')
                        if not os.path.exists(vrfile):
                            self.log.info("vasprun.xml file does not exist")
                        else:
                            if vac not in mater[func]:
                                mater[func][vac] = {}
                                
                            vr = Vasprun(vrfile)
                            gap,cbm,vbm,direct = vr.eigenvalue_band_properties
#                            mater[func][vac].update({"Etot": vr.final_energy})
#                            mater[func][vac].update({"Egap": gap})
                            mater[func][vac].update({"VBM": vbm})
                            
                            if vac == "vac_20":
                                ## get mu = energy per formula unit
                                ## usually vacuum spacing of 20 A is sufficiently well-converged
                                ## so we'll use that energy...
                                structure = Poscar.from_file(os.path.join(dir_vac,'POSCAR')).structure
                                formula_units = (structure.composition.num_atoms /
                                                 structure.composition.reduced_composition.num_atoms)
                                mater[func].update({"mu": vr.final_energy/formula_units})
                                mater[func].update({"Egap": gap})

                if "vac_20" not in myutils.listdironly(dir_func):
                    self.log.info("can't find the vac_20 subdirectory :( \
                                   I don't know what converged energy to use...")
                                    
                if eps_slab and d_slab:
                    mater[func].update({"eps_slab": eps_slab[i]})
                    mater[func].update({"d_slab": d_slab[i]}) 
               
            else:
                ## for bulk system, this should be a lot more straightforward
                structure = Poscar.from_file(os.path.join(dir_func,'POSCAR')).structure
                formula_units = (structure.composition.num_atoms /
                                 structure.composition.reduced_composition.num_atoms)
                
                vrfile = os.path.join(dir_func,'vasprun.xml')
                if not os.path.exists(vrfile):
                    self.log.info("vasprun.xml file does not exist")
                else:
                    vr = Vasprun(vrfile)
                    gap,cbm,vbm,direct = vr.eigenvalue_band_properties
                    mater[func].update({"mu": vr.final_energy/formula_units})
#                    mater[func].update({"Egap": gap})
#                    mater[func].update({"VBM": vbm})
            
         
        ## write json file    
        self.write_to_json(mater)
             
        return


    def create_entry_from_formula(self, mu_limit, formula):

        """
        create a database entry with the format 
        {
            functional: {
                mu_mu-limit (formula): xxx
            }
        } 
        the mu is calculated from a specified formula   
        
        """
        
        if self.dbentry:
            ## load existing json file
            mater = self.load_from_json()
        else:
            ## or initialize empty dict
            mater = {}
            self.dbentry = "%s.json"%self.system

        ## parse out the formula
        terms = []
        subtract = False
        for i,term in enumerate(re.split('(\+|-)',formula)):
            if term == "-":
                subtract = True
                continue
            elif subtract is True:
                terms.append((-float(term.split('*')[0]),term.split('*')[1]))
                subtract = False
            else:
                terms.append((float(term.split('*')[0]),term.split('*')[1]))
                
        for func in self.funcs:
          
            if func not in mater:
                mater[func] = {}
                
            mu = 0
            for (coeff,system) in terms:
                dbentry_system = os.path.join(self.dir_db,"%s.json"%system)
                if not os.path.exists(dbentry_system):
                    self.log.info("can't find the database entry for %s"%system)
                else:
                    with open(dbentry_system, 'r') as file:
                        system = json.loads(file.read())                
                    mu += coeff*system[func]["mu"]
            
            mater[func] = {"mu_%s (%s)"%(mu_limit, formula) : mu}
           
        ## write json file
        self.write_to_json(mater)
        
        return
        

    def write_to_json(self, mater):
        
        with open(os.path.join(self.dir_db,self.dbentry), 'w') as file:
             file.write(json.dumps(mater,sort_keys=True,indent=4))
                
        return

    
    def load_from_json(self):
        
        with open(os.path.join(self.dir_db, self.dbentry), 'r') as file:
            mater = json.loads(file.read())
                
        return mater
        

def main(args):
    
    ## define a main function callable from another python script
    parser = argparse.ArgumentParser(description='Create or modify simple database entries')
    parser.add_argument('main_system',help='the main system e.g. MoS2, WSe2')
    parser.add_argument('dir_db',help='path to the database directory')
    parser.add_argument('method',help='how to construct database entry \
                        ("from_vasp" or "from_formula")')
    
    parser.add_argument('--dir_dft',help='required if method is "from_vasp": \
                        path to main directory containing DFT calcs for this system')
    
    parser.add_argument('--mu_limit',help='required if method is "from_formula": \
                        which chemical potential limit to consider, e.g. Mo-rich')
    parser.add_argument('--formula',help='required if method is "formula": \
                        formula to use to calculate chemical potential, e.g. 0.5*MoS2-0.5*Mo')
    
    parser.add_argument('--monolayer',help='is this a monolayer system with vacuum dependence?',
                        default=False,action='store_true')
    parser.add_argument('--funcs', nargs='+', help='list each functional separated by a space, \
                        use same naming convention as subdirectories in dir_dft',
                        default = ["GGA", "SCAN_vdW"])
    parser.add_argument('--eps_slab', nargs='+', help='eps_slab calculated manually for the slab. \
                        List each dielectric constant separated by a space, \
                        in the order corresponding to funcs and d_slab.')
    parser.add_argument('--d_slab', nargs='+', help='d_slab calculated manually for the slab. \
                        List each slab thickness separated by a space, \
                        in the order corresponding to funcs and eps_slab.')
    parser.add_argument('--dbentry', help='existing database entry to append to')
    parser.add_argument('--logfile',help='logfile to save output to')
       
    ## read in the above arguments from command line
    args = parser.parse_args(args)
    
    ## set up logging
    if args.logfile:
        myLogger = myutils.setup_logging(args.logfile)
    else:
        myLogger = myutils.setup_logging()
     
        
    if args.method == "from_vasp":
        if not args.dir_dft:
            myLogger.info('For method "from_vasp", dir_dft is a required argument')
        else:
            db = DatabaseEntry(args.main_system, args.monolayer, args.dir_dft, 
                               args.dir_db, args.dbentry, args.funcs, myLogger)
            db.create_entry_from_vasp(eps_slab=args.eps_slab,
                                      d_slab=args.d_slab)
            
    elif args.method == "from_formula":
        if not args.mu_limit or not args.formula:
            myLogger.info('For method "from_formula", mu_limit and formula are required arguments')
        else:
            myLogger.info('For method "from_formula", mu_limit and formula are required arguments')
            db = DatabaseEntry(args.main_system, args.monolayer, args.dir_dft, 
                               args.dir_db, args.dbentry, args.funcs)     
            db.create_entry_from_formula(mu_limit=args.mu_limit, 
                                         formula=args.formula)
            
    else:
        myLogger.info('Invalid method!')

    
    
if __name__ == '__main__':
    

#    main(["--h"])    
#    main(["MoS2","minidb/","from_vasp","--dir_dft","/Users/annemarietan/ufrc/MoS2/monolayer_unitcell/eqm_a0/",
#          "--monolayer","--funcs","GGA","SCAN_vdW","--eps_slab","17.18","16.27","--d_slab","5.4","5.4"])
#    main(["Mo","minidb/","from_vasp","--dir_dft","/Users/annemarietan/orange/bccMo/",
#          "--funcs","SCAN_vdW","--dbentry","Mo.json"])
#    main(["S","minidb/","from_formula","--mu_limit","Mo-rich","--formula","0.5*MoS2-0.5*Mo"])
    
    
    main(sys.argv[1:]) 
    

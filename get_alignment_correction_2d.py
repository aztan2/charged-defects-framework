import os
import argparse
import logging
import time
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')


if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description='Estimate alignment correction.')
    parser.add_argument('vref',help='path to bulk LOCPOT file')
    parser.add_argument('vdef',help='path to defect LOCPOT file')
    parser.add_argument('encut',type=int,help='cutoff energy (eV)')
    parser.add_argument('q',type=int,help='charge (conventional units)')
    parser.add_argument('--vfile',help='vline .dat file',default='vline-eV.dat')
    parser.add_argument('--noplot',help='do not generate plots',default=False,action='store_true')
    parser.add_argument('--logfile',help='logfile to save output to')
       
    ## read in the above arguments from command line
    args = parser.parse_args()
    
    ## set up logging
    if args.logfile:
        logging.basicConfig(filename=args.logfile,filemode='w',
                            format='%(levelname)s:%(message)s',level=logging.DEBUG)    
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(logging.Formatter('%(levelname)s:%(message)s'))
        logging.getLogger('').addHandler(console)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.DEBUG)
    
    ## basic command to run sxdefectalign2d
    command = ['~/sxdefectalign2d', '--vasp',
               '--ecut', str(args.encut/13.6057), ## convert eV to Ry
               '--vref', args.vref,
               '--vdef', args.vdef]
    
    ## initialize the range of shift values bracketing the optimal shift
    smin, smax = -np.inf, np.inf
    shift = 0.0
    shifting = 'right'
    done = False
    counter = -1

    time0 = time.time()
    while not done and counter < 20:
        counter += 1
        ## run sxdefectalign2d with --shift <shift>
        command1 = command + ['--shift', str(shift), '--onlyProfile']
        os.system(' '.join(command1))
        
        ## read in the potential profiles from vline-eV.dat
        ## z  V^{model}  \DeltaV^{DFT}  V^{sr}
        data = np.loadtxt(args.vfile)
        
        ## assumes that the slab is in the center of the cell vertically!
        ## select datapoints corresponding to 2 bohrs at the top and bottom
        ## of the supercell (i.e. a total of 4 bohrs in the middle of vacuum)
        z1 = np.min([i for i,z in enumerate(data[:,0]) if z > 2.])
        z2 = np.min([i for i,z in enumerate(data[:,0]) if z > (data[-1,0]-2.)])
        
        ## fit straight lines through each subset of datapoints
        m1,C1 = np.polyfit(data[:z1,0],data[:z1,-1],1)
        m2,C2 = np.polyfit(data[z2:,0],data[z2:,-1],1)
        logging.debug("Slopes: %.8f %.8f; Intercepts: %.8f %.8f"%(m1,m2,C1,C2))
        
        ## check the slopes and intercepts of the lines
        ## and shift the charge along z until the lines are flat
        if abs(m1) < 1e-3 and abs(m2) < 1e-3 and abs(C1-C2) < 1e-3:
            done = True
            break
        elif m1*m2 < 0:
            logging.info("undetermined...make a tiny shift and try again")
            if shifting == 'right':
                shift += 0.01
            else: 
                shift -= 0.01
            logging.info("try shift = %.8f"%shift)
        elif (m1+m2)*np.sign(args.q) > 0:
            smin = shift
            if smax == np.inf:
                shift += 1.0
            else:
                shift = (smin+smax)/2.0
            shifting = 'right'
            logging.debug("optimal shift is in [%.8f, %.8f]"%(smin,smax))
            logging.info("shift charge in +z direction; try shift = %.8f"%shift)
        elif (m1+m2)*np.sign(args.q) < 0:
            smax = shift
            if smin == -np.inf:
                shift -= 1.0
            else:
                shift = (smin+smax)/2.0
            shifting = 'left'
            logging.debug("optimal shift is in [%.8f, %.8f]"%(smin,smax))
            logging.info("shift charge in -z direction; try shift = %.8f"%shift)
        
        if not args.noplot:
            plt.figure()
            plt.plot(data[:,0],data[:,2],'r',label=r'$V_{def}-V_{bulk}$')
            plt.plot(data[:,0],data[:,1],'g',label=r'$V_{model}$')
            plt.plot(data[:,0],data[:,-1],'b',label=r'$V_{def}-V_{bulk}-V_{model}$')
            plt.xlabel("distance along z axis (bohr)")
            plt.ylabel("potential (eV)")
            plt.xlim(data[0,0],data[-1,0])
            plt.legend() 
#            plt.show()
            plt.savefig(os.getcwd()+'/alignment.png')
                       
    if done:
        C_ave = (C1+C2)/2
        logging.info("DONE! shift = %.8f & alignment correction = %.8f"%(shift,C_ave))
        ## run sxdefectalign2d with --shift <shift> -C <C_ave> > correction
        command2 = command + ['--shift', str(shift),
                              '-C', str(C_ave),
                              '> correction']
        os.system(' '.join(command2))
    else:
        logging.info("Could not find optimal shift after 20 tries :(")
    
    logging.debug("Total time taken (s): %.2f"%(time.time()-time0))
    

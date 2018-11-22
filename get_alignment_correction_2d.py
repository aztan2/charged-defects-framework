import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
#plt.switch_backend('agg')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Estimate alignment correction.')
    parser.add_argument('q',type=int,help='charge (conventional units)')
    parser.add_argument('--vfile',help='vline .dat file',default='vline-eV.dat')
    parser.add_argument('--plot',help='whether to generate plots',default=True)
       
    ## read in the above arguments from command line
    args = parser.parse_args()
 
    data = np.loadtxt(args.vfile)
    
    ## assume that the vacuum is at the top/bottom of the supercell
    z1 = np.min([i for i,z in enumerate(data[:,0]) if z > 2.])
    z2 = np.min([i for i,z in enumerate(data[:,0]) if z > (data[-1,0]-2.)])
    m1,C1 = np.polyfit(data[:z1,0],data[:z1,-1],1)
    m2,C2 = np.polyfit(data[z2:,0],data[z2:,-1],1)
    print (m1,m2,C1,C2)
    if abs(m1) < 1e-3 and abs(m2) < 1e-3 and abs(C1-C2) < 1e-3:
        C_ave = (C1+C2)/2
        print ("done! alignment correction: %f"%C_ave)
    elif m1*m2 < 0:
        print ("undetermined...make a tiny shift and try again")
    elif abs(m1)*np.sign(args.q) > 0:
        print ("shift charge in +ve z direction")
    elif abs(m1)*np.sign(args.q) < 0:
        print ("shift charge in -ve z direction")
    
    if args.plot:
        plt.figure()
        plt.plot(data[:,0],data[:,2],'r',label=r'$V_{def}-V_{bulk}$')
        plt.plot(data[:,0],data[:,1],'g',label=r'$V_{model}$')
        plt.plot(data[:,0],data[:,-1],'b',label=r'$V_{def}-V_{bulk}-V^{LR}$')
        plt.xlabel("distance along z axis (bohr)")
        plt.ylabel("potential (eV)")
        plt.xlim(data[0,0],data[-1,0])
        plt.legend() 
#        plt.show()
        plt.savefig(os.getcwd()+'/alignment.png')
    
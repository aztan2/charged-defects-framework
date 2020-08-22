import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')


def get_skiprows(f):
    ## the data we are interested in starts after the '&' sign
    for linenum,line in enumerate(f):
        if line[0] == '&':
            return (linenum+1)

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Estimate alignment correction in bulk systems. \n'
                                     'You must already have run the sxdefectalign script once '
                                     'to generate the vline-eV-a0/1/2.dat files. \n'
                                     'Afterwards, run the sxdefectalign script again '
                                     'with the additional tag -C <averaged alignment correction>')
    parser.add_argument('defpos_a0',type=float,
                        help='position of defect along first axis (relative coords)')
    parser.add_argument('defpos_a1',type=float,
                        help='position of defect along second axis (relative coords)')
    parser.add_argument('defpos_a2',type=float,
                        help='position of defect along third axis (relative coords)')
    parser.add_argument('--noplot',help='do not generate plots',default=False,action='store_true')
       
    ## read in the above arguments from command line
    args = parser.parse_args()

    C_ave = 0
    for axis,defpos in zip(['a0','a1','a2'],[args.defpos_a0,args.defpos_a1,args.defpos_a2]):    
        with open('vline-eV-%s.dat'%axis,'r') as f:
            data = np.loadtxt('vline-eV-%s.dat'%axis,skiprows=get_skiprows(f))
        
        plat_pos = (defpos+0.5)%1 * data[-1,0]
        plat_v = [] 
        for x,v in zip(data[:,0],data[:,2]):
            if x > plat_pos - 1 and x < plat_pos + 1:
                plat_v.append(v)
        C = np.average(plat_v)
        C_ave += C
        print ("alignment correction: %f"%C)
    
        if not args.noplot: 
            plt.figure()
            plt.plot(data[:,0],data[:,2],'g',label=r'$V_{def}-V_{bulk}-V^{LR}$')
            plt.hlines(C,data[0,0],data[-1,0],'g',linestyle='dashed',
                       label="alignment correction: %f"%C)
            plt.xlabel("distance along %s axis (bohr)"%axis)
            plt.ylabel("potential (eV)")
            plt.xlim(data[0,0],data[-1,0])
            plt.legend() 
            plt.savefig(os.getcwd()+'/alignment_%s.png'%axis)
            
    print ("averaged alignment correction: %f"%(C_ave/3))
    
    
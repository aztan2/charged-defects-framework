import numpy as np
import pandas as pd
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar


if __name__ == '__main__':

    
    folder = "mp-2815_MoS2/monolayer_Svac/GGA/mag/"
    writer = pd.ExcelWriter('test.xlsx')
    
    dfs = []
    col_names =  ['vacuum', 'supercell', 'N', '1/N', 'E_def', 'E_bulk']

    qs = [-1,1]
    cells = [(3,3),(3,2),(4,4),(4,2),(5,5)]
    vacs = [10,15,20]
    
    ## set up dataframe for neutral defect first
    df0 = pd.DataFrame(columns = col_names)
    for vac in vacs:
        for cell in cells:            
            folder0 = folder + "charge_0/%dx%dx1/vac_%d/"%(cell[0],cell[1],vac)
            outcar = Outcar(folder0 + "OUTCAR")
            outcar_ref = Outcar(folder0 + "bulkref/OUTCAR")
            natoms = np.sum(Poscar.from_file(folder0 + "bulkref/POSCAR").natoms)
            df0.loc[len(df0)] = [vac,
                                 '%dx%dx1'%(cell[0],cell[1]),
                                 natoms,
                                 1/natoms,
                                 outcar.final_energy,
                                 outcar_ref.final_energy]
            ## add the parsing of chem pot, vbm
    df0.to_excel(writer,'charge_0')

    ## modify dataframe for charged defects
    for q in qs:
        df = df0.copy(deep=True)
        i = -1
        for vac in vacs:
            for cell in cells: 
                i += 1
                folder1 = folder + "charge_%d/%dx%dx1/vac_%d/"%(q,cell[0],cell[1],vac)
                outcar = Outcar(folder1 + "OUTCAR")
                df.loc[i,'E_def'] = outcar.final_energy
        df.to_excel(writer,'charge_%d'%q)
    
    writer.save()
                
#    df_1 = df_1.assign(new_col=df_1["E_def"] - df_1["E_bulk"])                
    
    
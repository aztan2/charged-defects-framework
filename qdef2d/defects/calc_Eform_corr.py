import os
import argparse
import pandas as pd


def calc(dir_def,xlfile):
    
    """
    Evaluate corrected defect formation energy.
    
    dir_def (str): path to the defect directory containing the excel file
    xlfile (str): excel filename to read/save the dataframe from/to
    
    """

    ## load list of dataframes from sheets from excel file    
    df = pd.read_excel(os.path.join(dir_def,xlfile),sheet_name=None)

   
    for q in [qi for qi in df.keys() if qi != 'charge_0']:           
        ## Finally, we can compute the corrected defect formation energy:
        ## Eform = Eform_uncorr + E_correction
        df[q]["E_form_corr"] = df[q].loc[:,'E_form_uncorr'] + df[q].loc[:,'E_corr']


    ## write the updated excel file
    writer = pd.ExcelWriter(os.path.join(dir_def,xlfile))
    for q in df.keys():  
        df[q].to_excel(writer, q, index=False)
    writer.save()   
    

if __name__ == '__main__':
    

    ## this script can also be run directly from the command line
    parser = argparse.ArgumentParser(description='Evaluate corrected defect formation energy.')
    parser.add_argument('dir_def',help='path to the defect directory containing the excel file')
    parser.add_argument('xlfile',help='excel filename to read/save the dataframe from/to')
       
    ## read in the above arguments from command line
    args = parser.parse_args()
    
    calc(args.dir_def, args.xlfile) 
    
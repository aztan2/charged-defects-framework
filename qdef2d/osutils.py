import os
import errno


def listdironly(path):
    
    return [d for d in os.listdir(path) 
            if os.path.isdir(os.path.join(path,d)) == True]


def joinpath(path,*paths):
    
    ## force the paths to have forward slash (unix-style)
    
    return os.path.join(path,*paths).replace('\\','/')


def check_file_exists(directory,filename):
    
    ## check if a file is present in a given directory
    
    files = [f for f in os.listdir(directory) 
                if f.startswith(filename)]
    if len(files) < 1:
        return False
    elif len(files) > 1:
        return False
    else:
        return True
   
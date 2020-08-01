import os
import logging


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
        print ("can't find %s file!"%filename)
        return False
    elif len(files) > 1:
        print ("more than 1 %s file found"%filename)
        return False
    else:
        return True
    
    
def setup_logging(logfile=None):

    logger = logging.getLogger()
    if logger.hasHandlers(): 
        logger.handlers = []
    logging.getLogger('matplotlib.font_manager').disabled = True
    
    if logfile:
        file = logging.FileHandler(logfile)
        file.setLevel(logging.DEBUG)
        file.setFormatter(logging.Formatter('%(levelname)s:%(message)s'))
        logger.addHandler(file) 
        
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(logging.Formatter('%(levelname)s:%(message)s'))
        logger.addHandler(console)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.DEBUG)
        
    return logger
        
        
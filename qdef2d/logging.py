import logging
    
    
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
        
        
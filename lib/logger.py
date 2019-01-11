# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 15:04:50 2018

@author: nicolas.chevrollier
"""
import os
import logging

def define_logger(path2output=None, outname=None, date=None):
    from logging.handlers import RotatingFileHandler
    
    logfilename = path2output+outname+'_'+date+'.log'
    print logfilename
    
    if not os.path.exists(logfilename):
        with open(logfilename, 'w') as f: pass
    
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG) 
    #formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
    formatter = logging.Formatter()
    file_handler = RotatingFileHandler(logfilename, 'a', 1000000, 1)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler) 
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    logger.addHandler(stream_handler)
    
    return logger

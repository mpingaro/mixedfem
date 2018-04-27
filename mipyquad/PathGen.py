# -*- coding: utf-8 -*-
"""
Created on 

@author: Dr. Marco Pingaro
@ref   : La Spienza University of Rome

Python Version 2.7
"""

import os

class PATH:
    def __init__(self):

        self.WS  = None           # Path Installation
        self.FOR = None           # Path Formulations
        self.FUN = None           # Path Functions
        self.IO = None            # Path Input/Output
        self.SOL = None           # Path Solver
        self.UT = None            # Path Utility

        self.CURRENT = None # Current path
        self.OUT = None  # Path Results


def GenPathFolders():
    # Call structure Path
    path = PATH()

    # Folders installation
    path.WS = "/home/mpingaro/temp/PyVEM" # Folder installation
    path.FOR = path.WS + "/formulations"
    path.FOR_SMALL = path.WS + "/formulations/small_displacement"
    path.FOR_LARGE = path.WS + "/formulations/large_displacement"
    path.FUN = path.WS + "/functions"
    path.IO = path.WS + "/io"
    path.SOL = path.WS + "/solver"
    path.UT = path.WS + "/utils"

    # Work Folders
    path.CURRENT = os.getcwd() # Current folder


    # Fix path Results, Figure and Geometry
    path.OUT = path.CURRENT + "/output/"

    if (os.path.exists(path.OUT)==False):
        # Gen folders
        os.system("mkdir " + path.CURRENT + "/output/")
    # TODO: not clear the results
    #else:

        # Clear folders
        #os.system("rm -rf "+ path.OUT + "/*")

    return path

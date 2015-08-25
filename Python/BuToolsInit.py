import os, sys

def BuToolsInit (verbose=False, checkInput=True, checkPrecision=1e-12):
#    btDir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    btDir = os.path.dirname(os.path.realpath(__file__))
    sys.path.append(btDir)
    print("Butools V2.0")
    packages = ['utils', 'mc', 'moments', 'reptrans', 'trace', 'ph', 'dph', 'map', 'dmap', 'fitting' , 'mam', 'queues']
    get_ipython().ex("import butools")
    get_ipython().ex("import numpy as np")
    get_ipython().ex("import numpy.matlib as ml")
    get_ipython().ex("import matplotlib.pyplot as plt")
    get_ipython().magic("%matplotlib inline")

    sys.stdout.write("Packages loaded: ")
    sys.stdout.flush()
    for p in range(len(packages)):
        get_ipython().ex("from butools.{0} import *".format(packages[p]))
        sys.stdout.write(packages[p])
        if p<len(packages)-1:
            sys.stdout.write(", ")
        else:
            sys.stdout.write("\n")
        sys.stdout.flush()
    
    if verbose==0 or verbose==False or verbose=="0" or verbose=="False":
        get_ipython().ex("butools.verbose = {0}".format(False))
    elif verbose==1 or verbose==True or verbose=="1" or verbose=="True":
        get_ipython().ex("butools.verbose = {0}".format(True))
    if checkInput==0 or checkInput==False or checkInput=="0" or checkInput=="False":      
        get_ipython().ex("butools.checkInput = {0}".format(False))
    elif checkInput==1 or checkInput==True or checkInput=="1" or checkInput=="True":
        get_ipython().ex("butools.checkInput = {0}".format(True))        
    get_ipython().ex("butools.checkPrecision = {0}".format(checkPrecision))        
    get_ipython().ex('print("Global variables: \\nbutools.verbose =", butools.verbose, ", butools.checkInput =", butools.checkInput, ", butools.checkPrecision =", butools.checkPrecision)')
    
if __name__=="__main__":
    if len(sys.argv)<2:
        BuToolsInit()
    elif len(sys.argv)<3:
        BuToolsInit(sys.argv[1])
    else:
        BuToolsInit(sys.argv[1], sys.argv[2])


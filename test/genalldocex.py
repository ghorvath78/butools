import sys
import os
from os.path import splitext,dirname,realpath

scriptDir = dirname(realpath(__file__))

for fName in os.listdir(scriptDir+'/gendocex'):
    if fName.endswith('.py'):
        os.system("python "+scriptDir+"/gendocex/"+fName)
    elif fName.endswith('.m'):
        os.system("cd {0}; matlab -nodesktop -nosplash -r \"addpath('/home/gabor/work/SMCSolver/QBD');addpath('/home/gabor/work/SMCSolver/MG1'); {1}; exit;\"; cd ..".format(scriptDir+"/gendocex",splitext(fName)[0]))
    elif fName.endswith('.ma'):
        os.system("math -script "+scriptDir+"/gendocex/"+fName)


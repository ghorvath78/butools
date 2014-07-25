import re
import os

def processMatlabDoc(fun, btModule, signature, doc):

    fname = "../Matlab/"+btModule+"/"+fun+".m";
    lines = []
    with open(fname, "rt") as f:
        lines = f.readlines()

    if len(lines)==0:
        print ("  File",fname,"not found!")

    newlines = ["%  " + signature + "\n", "%  \n"]
    for dl in doc:
        newlines.append("%  " + dl)
    newlines.append("\n")

    found = False    
    for l in range(len(lines)):
        lsplit = lines[l].split()
        if len(lsplit)>0 and lsplit[0]=="function":
            newlines.extend(lines[l:])
            found = True
            break
            
    if not found:
        print ("  Can not find the start of the function in the matlab code!")
            
    with open(fname, "wt") as f:
        f.writelines(newlines)

def processPythonDoc(fun, btModule, doc):
    
    dir = "../Python/butools/"+btModule;
    if not os.path.exists (dir):
        print("  Can not find the required python directory!")
        return
        
    pyfiles = [name for name in os.listdir(dir) if name.endswith('.py')]
    
    found = False    
    for py in pyfiles:
        fullName = dir+"/"+py
        with open(fullName, "rt") as f:
            lines = f.readlines()
            for l in range(len(lines)):
                lSplit = re.split(" |\t|\(", lines[l])
                if lSplit[0]=="def" and lSplit[1]==fun:
                    found = True
                    # auto-detect the indent
                    ind = ""
                    x = 1
                    while len(lines[l+x].strip())==0:
                        x+=1
                    for i in range(len(lines[l+x])):
                        if lines[l+x][i].isspace():
                            ind = ind + lines[l+x][i]
                        else:
                            break
                    # prepare docstring
                    docs = [ind + "\"\"\"\n"]
                    for dl in doc:
                        docs.append(ind + dl)
                    docs.append(ind + "\"\"\"\n")
                    # remove old docstring, if there is one
                    codestart = l+1
                    if lines[codestart].startswith(ind+"\"\"\""):
                        codestart+=1
                        while not lines[codestart].startswith(ind+"\"\"\""):
                            codestart+=1
                        codestart+=1
                    # new content for the file
                    nlines = lines[:l+1] + docs + lines[codestart:];
                    # write to file
                    with open(fullName, "wt") as o:
                        o.writelines(nlines)
                    break
        if found:
            break
    
    if not found:
        print ("  Can not find the start of the function in the python code!")


rstfiles = [os.path.splitext(name)[0] for name in os.listdir('source') if name.endswith('.rst') and name[0].isupper()]

for fun in rstfiles:
    print("processing ",fun+".rst", "...")
    with open("source/"+fun+".rst", "rt") as f:
        btModule = None
        hasMatlab = False
        hasPython = False
        # parse rst file
        lines = f.readlines()
        l = 0
        while l<len(lines):
            sLine = lines[l].strip()
            if sLine.startswith(".. currentmodule::"):
                btModule = sLine.split(" ")[2].split(".")[1]
            elif sLine.startswith("* - Matlab:"):
                hasMatlab = True
                signature = lines[l+1].split("`")[1]
            elif sLine.startswith("* - Python/Numpy:"):
                hasPython = True
            elif len(sLine)>0 and sLine[0].isalpha() and l>0:
                # copy all lines to doc, till the examples
                doc = []
                while l<len(lines) and lines[l].strip()!="Examples":
                    if len(lines[l])<4:
                        doc.append("\n")
                    else:
                        doc.append(lines[l][4:])
                    l+=1
                break;
            l+=1;
        # remove last empty lines
        while len(doc)>0 and doc[-1]=="\n":
            doc.pop()
        # remove :math:
        doc = [d.replace(":math:","") for d in doc]           
        # write documentation to Matlab files
        if hasMatlab:
            processMatlabDoc(fun, btModule, signature, doc)
        # write documentation to Python files
        if hasPython:
            processPythonDoc(fun, btModule, doc)
                
                

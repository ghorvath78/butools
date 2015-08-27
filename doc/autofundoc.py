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

def processMathematicaDoc(fun, btModule, mathematicaSignature):
    
    # open package file
    oneLineDescr = None
    with open("source/"+btModule+".rst", "rt") as f:
        lines = f.readlines()
        l = 0
        while l<len(lines):
            if re.match(".*:py:func:`{0} ".format(fun), lines[l]):
               oneLineDescr = re.match(".*- (.*)",lines[l+1]).groups()[0]
               break
            l+=1
    
    # find corresponding mathematica file and update the usage string
    if oneLineDescr!=None:    
        dir = "../Mathematica/BuTools"
        if not os.path.exists (dir):
            print("  Can not find the required mathematica directory!")
            return
        mfiles = [name for name in os.listdir(dir) if name.lower()==btModule+".m"]
        if len(mfiles)==0:
            print("  Can not find mathematica file for package "+btModule)
            return
        fullName = dir+"/"+mfiles[0]
        found = False    
        nlines=[]
        with open(fullName, "rt") as f:
            lines = f.readlines()
            for l in lines:
                if l.lstrip().startswith(fun+"::usage"):
                    nlines.append(fun+'::usage = "{0}: {1}";\n'.format(mathematicaSignature, oneLineDescr))
                    found = True
                else:
                    nlines.append(l)
        if not found:
            print("   Can not file usage line for function "+fun)
            return
        else:
            with open(fullName, "wt") as o:
                o.writelines(nlines)

def exampleForFunction(funName, output):
    docExFiles = [name for name in os.listdir('../test/docex') if name.endswith('_{0}.docex'.format(output))]
    strex = []    
    for fil in docExFiles:
        with open("../test/docex/"+fil, "rt") as f:
            lines = f.readlines()
            capture = False
            for l in lines:
                if capture and (l.startswith("=== ") or l.startswith('"=== ')):
                    capture = False
                elif capture:
                    strex.append(l)
                if l.startswith("=== {0} ===".format(funName)) or l.startswith('"=== {0} ==="'.format(funName)):
                    capture = True
    return "".join(strex)     

rstfiles = [os.path.splitext(name)[0] for name in os.listdir('rawsource') if name.endswith('.rst') and name[0].isupper()]
exStringHeader = [("matlab","For Matlab:"), ("mathematica","For Mathematica:"), ("python","For Python/Numpy:")]

for fun in rstfiles:
    print("processing ",fun+".rst", "...")
    newLines = []
    with open("rawsource/"+fun+".rst", "rt") as f:
        btModule = None
        hasMatlab = False
        hasPython = False
        hasMathematica = False
        # parse rst file
        lines = f.readlines()
        l = 0
        while l<len(lines):
            sLine = lines[l].strip()
            if sLine.startswith(".. currentmodule::"):
                btModule = sLine.split(" ")[2].split(".")[1]
            elif sLine.startswith("* - Matlab:"):
                hasMatlab = True
                matlabSignature = lines[l+1].split("`")[1]
            elif sLine.startswith("* - Python/Numpy:"):
                hasPython = True
            elif sLine.startswith("* - Mathematica:"):
                hasMathematica = True
                mathematicaSignature = lines[l+1].split("`")[1]
            elif len(sLine)>0 and sLine[0].isalpha() and l>0:
                # copy all lines to doc, till the examples
                doc = []
                while l<len(lines) and lines[l].strip()!="Examples":
                    if len(lines[l])<4:
                        doc.append("\n")
                    else:
                        doc.append(lines[l][4:])
                    newLines.append(lines[l])
                    l+=1
                break
            newLines.append(lines[l])
            l+=1
        # remove last empty lines
        while len(doc)>0 and doc[-1]=="\n":
            doc.pop()
        # remove :math:
        doc = [d.replace(":math:","") for d in doc]           
        # write documentation to Matlab files
        if hasMatlab:
            processMatlabDoc(fun, btModule, matlabSignature, doc)
        # write documentation to Python files
        if hasPython:
            processPythonDoc(fun, btModule, doc)
        # write documentation to Mathematica files
        if hasMathematica:
            processMathematicaDoc(fun, btModule, mathematicaSignature)
                
        # append examples
        exString = [exampleForFunction(fun, es[0]) for es in exStringHeader]
        if any(map(len,exString)):
            newLines.append("    Examples\n")
            newLines.append("    ========\n")
            for i in range(len(exString)):
                newLines.append("    "+exStringHeader[i][1]+"\n")
                newLines.append("\n")
                newLines.append("    "+"    ".join(exString[i].splitlines(True)))
                newLines.append("\n")
        else:
            newLines = lines
        with open("source/"+fun+".rst", "wt") as o:
            o.writelines(newLines)
        

# -*- coding: utf-8 -*-
"""
Created on Tue May 19 14:15:46 2015

@author: gabor
"""

import io
import sys
import re
from contextlib import redirect_stdout

# open file

def processTestInit(funName):
    if forDocEx:
        if output=="MATLAB":
            print("disp('=== "+funName+" ===');")
        elif output=="Mathematica":
            print('Print["=== '+funName+' ==="//OutputForm];')
        elif output=="Python":
            print("print('=== "+funName+" ===')")
    else:
        if output=="MATLAB":
            print("disp('{0}')".format("="*40))
            print("help {0}".format(funName))
        elif output=="Mathematica":
            print('Print["{0}"];'.format("="*40))
            print('?"{0}";'.format(funName))
        elif output=="Python":
            print('print("{0}")'.format("="*40))
            print("print('Testing BuTools function {0}')".format(funName))    

def printQuotedMsg(what):
    if output=="MATLAB":
        what = re.compile("'").sub("''", what)
        print("disp('{0}');".format(what))
    elif output=="Mathematica":
        what = re.compile('"').sub('\\"', what)
        print('Print["{0}"//OutputForm];'.format(what))
    elif output=="Python":
        what = re.compile("'").sub('"', what)
        print("print('{0}')".format(what))

def printMsg(what):
    if output=="MATLAB":
        text = "disp({0});".format(what)
    elif output=="Mathematica":
        text = 'Print[{0}];'.format(what)
    elif output=="Python":
        text = "print({0})".format(what)
    if forDocEx:
        printQuotedMsg(">>> " + text)
    print(text)
        
def printVar (varName, echo):
    if echo:
        if not forDocEx:
            printQuotedMsg("{0} = ".format(varName))
        printMsg(varName)

def formatVec(varName, entries):
    if output=="MATLAB":
        return varName + " = [" + ", ".join(entries) + "];"
    elif output=="Mathematica":
        return varName + " = {" + ", ".join(entries) + "};"
    elif output=="Python":
        return varName + " = ml.matrix([[" + ", ".join(entries) + "]])"

def formatList(varName, entries):
    if output=="MATLAB":
        return varName + " = [" + ", ".join(entries) + "];"
    elif output=="Mathematica":
        return varName + " = {" + ", ".join(entries) + "};"
    elif output=="Python":
        return varName + " = [" + ", ".join(entries) + "]"

def formatMat(varName, entries):
    if output=="MATLAB":
        return varName + " = [" + "; ".join([", ".join(x) for x in entries]) + "];"
    elif output=="Mathematica":
        return varName + " = {{" + "},{".join([", ".join(x) for x in entries]) + "}};"
    elif output=="Python":
        return varName + " = ml.matrix([[" + "],[".join([", ".join(x) for x in entries]) + "]])"            

def formatRange(varName, rfrom, rto, rstep):
    if output=="MATLAB":
        return "{0} = ({1}:{2}:{3});".format(varName, rfrom, rstep, rto)
    elif output=="Mathematica":
        return "{0} = Range[{1},{2},{3}];".format(varName, rfrom, rto, rstep)
    elif output=="Python":
        return "{0} = np.arange({1},{2},{3})".format(varName, rfrom, str(float(rto)+float(rstep)), rstep)

def formatVar(varName, value):
    if output=="MATLAB":
        return varName + " = " + value + ";"
    elif output=="Mathematica":
        return varName + " = " + value + ";"
    elif output=="Python":
        return varName + " = " + value
        
def formatPlot(x,y):
    if output=="MATLAB":
        return "plot({0},{1});".format(x,y)
    elif output=="Mathematica":
        return "ListLinePlot[Transpose[{{{0}, {1}}}]]".format(x,y)
    elif output=="Python":
        return "plt.plot({0},{1})".format(x,y)    

def printAssert(cond, errMsg):
    if output=="MATLAB":
        print("assert({0}, '{1}');".format(cond, errMsg))
    elif output=="Mathematica":
        print('Assert[{0}, "{1}"];'.format(cond, errMsg))
    elif output=="Python":
        print("assert {0}, '{1}'".format(cond, errMsg))

def convertCode(code):
    if output=="Mathematica":
        # fix parenthesis
        code = re.compile('NPArray\[([^\]]+)\]').sub('\\1', code)
        code = re.compile('\./').sub('/', code)
        return code
    elif output=="MATLAB":
        # replace few function names
        code = re.compile('NPArray\[([^\]]+)\]').sub('\\1', code)
        code = re.compile('True').sub('true', code)
        code = re.compile('False').sub('false', code)
        code = re.compile('([a-zA-Z])([.])([a-zA-Z])').sub('\\1*\\3', code)
        code = re.compile('Norm\[').sub('norm(', code)
        code = re.compile('Abs\[').sub('abs(', code)
        code = re.compile('Max\[').sub('max(', code)
        code = re.compile('Min\[').sub('min(', code)
        code = re.compile('All\[').sub('all(', code)
        code = re.compile('Diff\[').sub('diff(', code)
        code = re.compile('Total\[').sub('sum(', code)
        code = re.compile('Length\[').sub('length(', code)
        code = re.compile('Infinity').sub('inf', code)
        # fix array indices
        code = re.compile('\[\[1\]\]').sub('', code)
        code = re.compile('\[\[([\w-]+)\]\]').sub('(\\1)', code)
        # fix parenthesis
        code = re.compile('\[').sub('(', code)
        code = re.compile('\]').sub(')', code)       
        code = re.compile('"').sub("'", code)
        return code
    elif output=="Python":
        # replace few function names
        code = re.compile('\./').sub('/', code)
        code = re.compile('([a-zA-Z])([.])([a-zA-Z])').sub('\\1*\\3', code)
        code = re.compile('NPArray\[([^\]]+)\]').sub('np.array(\\1)', code)
        code = re.compile('Norm\[').sub('la.norm(', code)
        code = re.compile('Abs\[').sub('np.abs(', code)
        code = re.compile('Max\[').sub('np.max(', code)
        code = re.compile('Min\[').sub('np.min(', code)
        code = re.compile('All\[').sub('np.all(', code)
        code = re.compile('Diff\[').sub('np.diff(', code)
        code = re.compile('Dot\[').sub('np.dot(', code)
        code = re.compile('Total\[').sub('np.sum(', code)
        code = re.compile('Length\[').sub('Length(', code)
        code = re.compile('Infinity').sub('np.inf', code)
        code = re.compile('&&').sub(' and ', code)
        code = re.compile('\|\|').sub(' or ', code)
        code = re.compile('\^').sub('**', code)
        # fix parenthesis
        code = re.compile('\[').sub('(', code)
        code = re.compile('\]').sub(')', code)
        # fix array indices
        code = re.compile('\(\(1\)\)').sub('[0]', code)
        code = re.compile('\(\(2\)\)').sub('[1]', code)
        return code
    

if sys.argv[1].startswith("MATLAB"):
    output = "MATLAB"
    funRetLeft, funRetRight = "[", "]"
elif sys.argv[1].startswith("Mathematica"):
    output = "Mathematica"
    funRetLeft, funRetRight = "{", "}"
elif sys.argv[1].startswith("Python"):
    output = "Python"
    funRetLeft, funRetRight = "", ""

forDocEx = sys.argv[1].endswith("DocEx")

usedVars = []
fio = io.StringIO()
with redirect_stdout(fio):
    with open(sys.argv[2],'r') as f:
        for line in f:
            line = line.lstrip()
            if not line or line.startswith("//"): 
                continue
            if line.startswith("package"):
                packName = line.split()[1]
                printQuotedMsg("---BuTools: {0} package test file---".format(packName))
                printQuotedMsg("Enable the verbose messages with the BuToolsVerbose flag")
                if output=="MATLAB":
                    print("global BuToolsVerbose;")
                    print("BuToolsVerbose = true;")
                elif output=="Mathematica":
                    print("BuTools`Verbose = True;")
                elif output=="Python":
                    print("butools.verbose = True")
                printQuotedMsg("Enable input parameter checking with the BuToolsCheckInput flag")
                if output=="MATLAB":
                    print("global BuToolsCheckInput;")
                    print("BuToolsCheckInput = true;")
                    if forDocEx:
                        print("format compact")
                        print("diary('{0}_matlab.docex');".format(packName))
                elif output=="Mathematica":
                    print("BuTools`CheckInput = True;")
                    print("On[Assert];")
                    if forDocEx:
                        print("tmpOut = $Output;")
                        print('stream = OpenWrite["{0}_mathematica.docex", FormatType -> InputForm, PageWidth -> Infinity];'.format(packName))
                        print('$Output = {stream};')
                        usedVars.append("stream")
                        usedVars.append("tmpOut")
                elif output=="Python":
                    print("butools.checkInput = True")
                    print("np.set_printoptions(precision=5,linewidth=1024)")
                    if forDocEx:
                        print("outFile = open('{0}_python.docex','w')".format(packName))
                        print("with redirect_stdout(outFile):")                        
            elif line.startswith("inputmsg"):
                if forDocEx:
                    continue;
                printQuotedMsg("Input:")
                printQuotedMsg("-" * 6)       
            elif line.startswith("testmsg"):
                if forDocEx:
                    continue;
                printQuotedMsg("Test:")
                printQuotedMsg("-" * 5)       
            elif line.startswith("test"):
                funName = line.split()[1]
                processTestInit(funName)
            elif line.startswith("defvec"):
                echo = line[6]=="@"
                varName, vectorString = re.match(r'\w+@?\s+(\w+)\s+#\s+(.*)',line).groups()
                entries = re.findall(r'[\s,]*([\w\.\-/]+)[\s,]*', vectorString)
                line = formatVec(varName, entries)
                if forDocEx:                
                    printQuotedMsg(">>> " + line)
                print(line)
                printVar(varName, echo and not forDocEx)
                usedVars.append(varName)
            elif line.startswith("deflist"):
                echo = line[7]=="@"
                varName, vectorString = re.match(r'\w+@?\s+(\w+)\s+#\s+(.*)',line).groups()
                entries = re.findall(r'[\s,]*([\w\.\-/]+)[\s,]*', vectorString)
                line = formatList(varName, entries)
                if forDocEx:                
                    printQuotedMsg(">>> " + line)
                print(line)
                printVar(varName, echo and not forDocEx)       
                usedVars.append(varName)
            elif line.startswith("defmat"):
                echo = line[6]=="@"
                varName, rowsString = re.match(r'\w+@?\s+(\w+)\s+#\s+\{(.*)\}',line).groups()
                rowsList = re.findall(r'\{([^\}]*)\}', rowsString)
                entries = [re.findall(r'[\s,]*([\w\.\-]+)[\s,]*', x) for x in rowsList]
                line = formatMat(varName, entries)
                if forDocEx:                
                    printQuotedMsg(">>> " + line)
                print(line)
                printVar(varName, echo and not forDocEx)       
                usedVars.append(varName)
            elif line.startswith("defrange"):
                echo = line[8]=="@"
                varName, rfrom, rto, rstep = re.match(r'\w+@?\s+(\w+)\s*#\s*([\w\.\-]+)\s*#\s*([\w\.\-]+)\s*#\s*([\w\.\-]+)\s*',line).groups()
                line = formatRange(varName, rfrom, rto, rstep)
                if forDocEx:                
                    printQuotedMsg(">>> " + line)
                print(line)
                printVar(varName, echo and not forDocEx)       
                usedVars.append(varName)
            elif line.startswith("defvar"):
                echo = line[6]=="@"
                varName, valueString = re.match(r'\w+@?\s+(\w+)\s+#\s+(.*)',line).groups()
                line = formatVar(varName, valueString)
                if forDocEx:                
                    printQuotedMsg(">>> " + line)
                print(line)
                printVar(varName, echo and not forDocEx)
                usedVars.append(varName)
            elif line.startswith("plot"):
                if not forDocEx:
                    continue
                x, y = re.match(r'plot\s+([\w\.\-]+)\s*#\s*([\w\.\-]+)\s*',line).groups()
                line = formatPlot(x,y)
                printQuotedMsg(">>> " + line)
            elif line.startswith("code"):                           
                echo = (line[4]=="@" or line[5]=="@")
                intro = (line[4]=="!" or line[5]=="!")
                firstnonws = 0
                while not line[firstnonws].isspace():
                    firstnonws+=1                    
                code = line[firstnonws:-1].strip()
                ccode = convertCode(code)
                outV = re.match(r'\{?([^\}=]+)\}?\s*=\s*(.*)',ccode)
                varLst = []
                if outV!=None:
                    varLst = re.findall(r'[\s,]*([\w\.\-]+)[\s,]*', outV.groups()[0])
                    usedVars += varLst
                    if len(varLst)==1:
                        ccode = "{0} = {1}".format(varLst[0], outV.groups()[1])
                    else:
                        ccode = "{0}{1}{2} = {3}".format(funRetLeft, outV.groups()[0], funRetRight, outV.groups()[1])
                if intro and not forDocEx:
                    printQuotedMsg(ccode+":")
                if output!="Python":
                    ccode += ";"
                if forDocEx:                
                    printQuotedMsg(">>> " + ccode)
                print(ccode)
                if echo:
                    for oa in varLst:
                        printVar(oa, echo)        
            elif line.startswith("print"):
                intro = False
                rest = line[6:-1]
                if line[5]=="!":
                    intro = True
                    rest = line[7:-1]
                strs = [x.strip() for x in rest.split("#")]
                for stri in strs:
                    if stri[0]!='"':
                        code = convertCode(stri)
                        if intro and not forDocEx:
                            printQuotedMsg(code+":")
                        printMsg(code)
                    elif not forDocEx:
                        printQuotedMsg(stri[1:-1])
            elif line.startswith("assert"):
                if forDocEx:
                    continue;
                forMathematica, forMatlab, forPython = False, False, False                
                for i in range(6,9):
                    if line[i]==' ':
                        break
                    elif line[i]=='W':
                        forMathematica = True
                    elif line[i]=='M':
                        forMatlab = True
                    elif line[i]=='P':
                        forPython = True
                if forMathematica==False and  forMatlab==False and forPython==False:
                    forMathematica, forMatlab, forPython = True, True, True
                if (output=="Mathematica" and forMathematica) or (output=="MATLAB" and forMatlab) or (output=="Python" and forPython):
                    cond, errMsg = re.match(r'\w+\s+([^#]+)#\s+(.*)',line).groups()
                    ccond = convertCode(cond)
                    printAssert(ccond, errMsg)
            else:
                raise Exception("unrecognized line "+line)

srcFile = fio.getvalue()

if forDocEx:
    if output=="MATLAB":
        srcFile += "diary off;\n"
        funName = "{0}MatlabGendocex".format(packName)
    elif output=="Python":
        # add indent
        lines = srcFile.splitlines(True)
        newLines = []
        state=0
        for l in lines:
            if state==1:
                newLines.append("\t" + l)
            else:
                newLines.append(l)
            if l.startswith("with redirect_stdout"):
                state=1                            
        newLines.append("outFile.close()\n")
        srcFile = "".join(newLines)        
        funName = "{0}PythonGendocex".format(packName)
    elif output=="Mathematica":
        srcFile = 'Unprotect[Print];\nPrint[args___] := Block[{$inMsg = True, result, str},\n   If[MatrixQ[args],\n    str = "{";\n    Do[str = StringJoin[str, ToString[args[[r]], FormatType -> InputForm]]; \n     If[r < Length[args], str = StringJoin[str, ",\\n "]], {r, Length[args]}];\n    str = StringJoin[str, "}"]; Print[str//OutputForm],\n    result = Print[args],\n    result = Print[args]\n    ]] /; ! TrueQ[$inMsg];\n'+srcFile
        srcFile += '\n$Output = tmpOut;\nClose[stream];\n'
        funName = "{0}MathematicaGendocex".format(packName)
else:
    funName = "Test{0}Package".format(packName)


if output=="Mathematica":
    print("{0}[]:=\nModule[{{{1}}},\n\t{2}\n];".format(funName, ", ".join(set(usedVars)), '\t'.join(srcFile.splitlines(True))))
elif output=="MATLAB":
    print("function {0} ()\n\t{1}end".format(funName, '\t'.join(srcFile.splitlines(True))))
elif output=="Python":
    print('import numpy as np\nimport numpy.matlib as ml\nimport matplotlib.pyplot as plt\nimport butools\nfrom butools.utils import *\nfrom butools.ph import *\nfrom butools.moments import *\nfrom butools.mc import *\nfrom contextlib import redirect_stdout\n')
    print("def {0}():\n\t{1}".format(funName, '\t'.join(srcFile.splitlines(True))))
    print('if __name__ == "__main__":\n\t{0}()'.format(funName))
    
    
                    
    
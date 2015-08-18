from antlr4 import *
from bttestLexer import bttestLexer
from bttestParser import bttestParser
from bttestListener import bttestListener
import sys 
import re
import os
from os.path import basename,splitext,dirname,realpath
from contextlib import redirect_stdout

class TestWriterBase(bttestListener):
    
    evals = {}
    nums = {}
    varLists = {}
    indent = 0
    forDepth = 0
    forBuffer = ""
    allForEnds = False
    docExEnds=False


    def printLine(self, what):
        if not forDocEx or not self.docExEnds:
            print(' '*self.indent + what)

    def printVar (self, enable, varName):
        if enable:
            if not forDocEx:
                self.printQuotedMsg("{0} = ".format(varName))
            self.printMsg(varName)
        
    def enterTestcase (self, ctx):
        self.docExEnds = False
        if forDocEx:
            print(' '*self.indent + self.strPrintCmd.format("=== "+str(ctx.ID())+" ==="))
        else:
            print(' '*self.indent + self.strPrintCmd.format("="*40))
            print(' '*self.indent + self.strPrintCmd.format("Testing BuTools function {0}").format(ctx.ID()))

    def enterTestmsg (self, ctx):
        if not forDocEx:
            self.printQuotedMsg("Test:")
            self.printQuotedMsg("-" * 5)       
    
    def enterInputmsg (self, ctx):
        if not forDocEx:
            self.printQuotedMsg("Input:")
            self.printQuotedMsg("-" * 6)       

    def enterEnddocex (self, ctx):
        self.docExEnds = True

    def exitDefvar (self, ctx):
        varName = str(ctx.ID())
        line = varName+" = " + self.evals[ctx.expression()] + self.lineDelim
        if forDocEx:                
            self.printQuotedMsg(">>> " + line)
        self.printLine(line)
        self.printVar (ctx.WBRESULT()!=None and not forDocEx, varName)

    def exitDefvec (self, ctx):
        varName = str(ctx.ID())
        numList = [self.evals[e] for e in ctx.vecspec().expression()]
        line = varName+" = "+self.vectorPre + ",".join(numList) + self.vectorPost + self.lineDelim
        if forDocEx:                
            self.printQuotedMsg(">>> " + line)
        self.printLine(line)
        self.printVar (ctx.WBRESULT()!=None and not forDocEx, varName)

    def exitDeflist (self, ctx):
        varName = str(ctx.ID())
        numList = [self.evals[e] for e in ctx.vecspec().expression()]
        line = varName+" = " + self.listPre + ", ".join(numList) + self.listPost + self.lineDelim
        if forDocEx:                
            self.printQuotedMsg(">>> " + line)
        self.printLine(line)
        self.printVar (ctx.WBRESULT()!=None and not forDocEx, varName)

    def exitDefmat (self, ctx):
        varName = str(ctx.ID())
        vecLists = ctx.vecspec()
        vecStrs = [", ".join([self.evals[e] for e in v.expression()])  for v in vecLists]
        line = varName+" = " + self.matrixPre + self.matrixColBreak.join(vecStrs) + self.matrixPost + self.lineDelim
        if forDocEx:                
            self.printQuotedMsg(">>> " + line)
        self.printLine(line)
        self.printVar (ctx.WBRESULT()!=None and not forDocEx, varName)

    def exitDefrange (self, ctx):
        varName = str(ctx.ID())
        rFrom, rTo, rStep = map(float,map(str,ctx.NUMBER()))
        line = varName + " = " + self.rangeCmd.format(rFrom, rTo, rStep) + self.lineDelim
        if forDocEx:                
            self.printQuotedMsg(">>> " + line)
        self.printLine(line)
        self.printVar (ctx.WBRESULT()!=None and not forDocEx, varName)
        
    def exitTprint (self, ctx):
        announce = str(ctx.ANNOUNCE())=="!"
        for e in ctx.expression():
            s = self.evals[e]
            if len(s)>0 and s[0]!='"' and announce and not forDocEx:
                self.printQuotedMsg(s+":")
            if s[0]!='"' or not forDocEx:
                self.printMsg(s)

    def exitTassert (self, ctx):
        if not forDocEx:
            print(' '*self.indent + self.assertCmd.format(self.evals[ctx.expression()], ctx.STRING().getText()[1:-1]))

    def exitArrindex (self, ctx):
        if ctx.NUMBER()!=None:
            self.nums[ctx] = int(str(ctx.NUMBER()))
        else:
            self.evals[ctx] = self.evals[ctx.expression()]                
        
    def exitArgumentList (self, ctx):
        self.evals[ctx] = ", ".join([self.evals[e] for e in ctx.expression()])

    def exitVarlist (self, ctx):
        self.varLists[ctx] = [self.evals[e] for e in ctx.expression()]

    def exitTlist (self, ctx):
        self.varLists[ctx] = [self.evals[e] for e in ctx.expression()]
        self.evals[ctx] = self.cellPre + ", ".join([self.evals[e] for e in ctx.expression()]) + self.cellPost

    def enterTplot (self, ctx):
        line = self.plotCmd.format(", ".join(map(str,ctx.ID())))
        if forDocEx and not self.docExEnds:
            self.printQuotedMsg(">>> " + line)
        else:
            self.printLine(line)
    
    def exitTcode (self, ctx):
        codeStr = ""
        varList = []
        toList = False        
        if ctx.varlist()!=None:
            if ctx.varlist().br==None and ctx.varlist().expression(0).tlist()!=None:
                toList = True
                varList = self.varLists[ctx.varlist().expression(0).tlist()]
            else:
                varList = self.varLists[ctx.varlist()]
        for i in range(len(varList)):
            if varList[i] in self.funMap:
                varList[i] = self.funMap[varList[i]]
        if len(varList)==1:
            codeStr += varList[0] + " = "
        elif len(varList)>1:
            codeStr += self.unPackOutPre + ", ".join(varList) + self.unPackOutPost + " = "
        codeStr += self.evals[ctx.expression()]
        if len(varList)>1 and toList:
            codeStr += self.listExtract
        if ctx.expression().ID()==None or (ctx.expression().ID().getText()!="For" and ctx.expression().ID().getText()!="EndFor"):
            codeStr += self.lineDelim
        if len(codeStr)>0:
            if self.forDepth>0:
                self.forBuffer += '\n' + ' '*(4*self.forDepth) + codeStr
            else:            
                if ctx.ANNOUNCE()!=None and not forDocEx and ctx.SILENT()==None:
                    for ln in codeStr.splitlines():
                        self.printQuotedMsg (ln + ":")
                if forDocEx and ctx.SILENT()==None:                
                    for ln in codeStr.splitlines():
                        self.printQuotedMsg(">>> " + ln)
                for ln in codeStr.splitlines():
                    self.printLine(ln)
            for v in varList:
                self.printVar (ctx.WBRESULT()!=None, v)                

# ==========================================================================================
class PythonWriter(TestWriterBase):

    funMap = {"Norm":"la.norm", "Abs":"np.abs", "Max":"np.max", "Min":"np.min", "Total":"np.sum",\
              "Diff":"np.diff", "ToArray":"np.array", "Eigenvalues":"la.eigvals", \
              "BuTools`CheckPrecision":"butools.checkPrecision", "Mean":"np.mean", "Inv":"la.inv",\
              "ExpMat":"la.expm","Sum":"np.sum","Ceil":"math.ceil","Linspace":"np.linspace",\
              "Eye":"ml.eye","BuTools`Verbose":"butools.verbose",\
              "Kron":"ml.kron","Diag":"Diag","Pinv":"la.pinv","MatrixMin":"np.min",\
              "MatrixMax":"np.max"}
    before = """\
import math
import numpy as np
import numpy.matlib as ml
import matplotlib.pyplot as plt
import butools
from butools.utils import *
from butools.ph import *
from butools.dph import *
from butools.map import *
from butools.moments import *
from butools.reptrans import*
from butools.mc import *
from butools.dmap import *
from butools.trace import *
from butools.mam import *
from butools.queues import *
from butools.fitting import *
from contextlib import redirect_stdout
import os

"""
    after = ""
    indent = 0
    printCmd = "print({0})"
    strPrintCmd = "print('{0}')"
    stringDelim = '"'
    vectorPre = "ml.matrix([["
    vectorPost = "]])"
    matrixPre = "ml.matrix([["
    matrixPost = "]])"
    matrixColBreak = "],["
    listPre = "["
    listPost = "]"
    cellPre = "["
    cellPost = "]"
    rangeCmd = "np.arange({0},{1},{2})"
    assertCmd = 'assert {0}, "{1}"'
    plotCmd =   "plt.plot({0})"
    listExtract = ""
    lineDelim = ""
    unPackOutPre = ""
    unPackOutPost = ""

    def printQuotedMsg(self, what):
        what = re.compile("'").sub('"', what)
        self.printLine("print('{0}')".format(what))

    def printMsg(self, what):
        text = "print({0})".format(what)
        if forDocEx:
            self.printQuotedMsg(">>> " + text)
        self.printLine(text)

    def enterTestfile (self, ctx):
        self.packName = str(ctx.ID())
        print('import sys')
        print('sys.path.append("'+os.path.realpath(os.path.dirname(os.path.realpath(__file__))+'/../Python")'))
        docExDir = os.path.realpath(os.path.dirname(os.path.realpath(__file__))+'/docex')+"/"
        print(self.before.format(self.packName))
        self.printQuotedMsg("---BuTools: {0} package test file---".format(self.packName))
        self.printQuotedMsg("Enable the verbose messages with the BuToolsVerbose flag")
        print(' '*self.indent + "butools.verbose = True")
        self.printQuotedMsg("Enable input parameter checking with the BuToolsCheckInput flag")
        print(' '*self.indent + "butools.checkInput = True")
        print(' '*self.indent + "np.set_printoptions(precision=5,linewidth=1024)")
        if forDocEx:
            print(' '*self.indent + "outFile = open('{0}_python.docex','w')".format(docExDir+self.packName))
            print(' '*self.indent + "with redirect_stdout(outFile):")
            self.indent += 4

    def exitTestfile (self, ctx):
        print(self.after.format(self.packName))

    def exitDefrange (self, ctx):
        varName = str(ctx.ID())
        rFrom, rTo, rStep = map(float,map(str,ctx.NUMBER()))
        line = varName + " = " + self.rangeCmd.format(rFrom, str(float(rTo)+float(rStep)), rStep) + self.lineDelim
        if forDocEx:                
            self.printQuotedMsg(">>> " + line)
        self.printLine(line)
        self.printVar (ctx.WBRESULT()!=None and not forDocEx, varName)

    def exitExpression (self, ctx):
        if ctx.primary()!=None:
            self.evals[ctx] = self.evals[ctx.primary()]            
        elif ctx.biop!=None:
            oper = ctx.biop.text
            if oper=="." or oper==".*":
                oper="*"
            elif oper=="^" or oper==".^":
                oper="**"
            elif oper=="./":
                oper="/"
            elif oper=="&&":
                oper=" and "
            elif oper=="||":
                oper=" or "
            self.evals[ctx] = self.evals[ctx.expression(0)]+oper+ self.evals[ctx.expression(1)]
        elif ctx.unop!=None:
            self.evals[ctx] = ctx.unop.text + self.evals[ctx.expression(0)]
        elif ctx.par1!=None:
            self.evals[ctx] = "(" + self.evals[ctx.expression(0)]+")"
        elif ctx.arrix!=None or ctx.lstix!=None:
            self.evals[ctx] = self.evals[ctx.expression(0)] + "[" + ", ".join([self.evals[e] for e in ctx.arrindexexpr()]) + "]"
        elif ctx.tlist()!=None:
            self.evals[ctx] = self.evals[ctx.tlist()]
        elif ctx.funpar!=None:
            # function call!!!
            funName = ctx.ID().getText()
            if funName=="AllNonNegative":
                self.evals[ctx] = "np.all({0}>=0)".format(self.evals[ctx.argumentList()])
            elif funName=="AllNegative":
                self.evals[ctx] = "np.all({0}<0)".format(self.evals[ctx.argumentList()])
            elif funName=="AllGreaterThanm1":
                self.evals[ctx] = "np.all({0}>=-1)".format(self.evals[ctx.argumentList()])
            elif funName=="AllLessThan1":
                self.evals[ctx] = "np.all({0}<=1)".format(self.evals[ctx.argumentList()])
            elif funName=="AllPositiveMatrix":
                self.evals[ctx] = "np.all({0}>0)".format(self.evals[ctx.argumentList()])
            elif funName=="AllNonNegativeMatrix":
                self.evals[ctx] = "np.all({0}>=0)".format(self.evals[ctx.argumentList()])
            elif funName=="AllNotGreaterThan1Matrix":
                self.evals[ctx] = "np.all({0}<=1)".format(self.evals[ctx.argumentList()])
            elif funName=="Dot":
                self.evals[ctx] = "{0}.dot({1})".format(self.evals[ctx.argumentList().expression(0)],self.evals[ctx.argumentList().expression(1)])
            elif funName=="Dim1":
                self.evals[ctx] = "{0}.shape[0]".format(self.evals[ctx.argumentList().expression(0)])
            elif funName=="RowSum":
                self.evals[ctx] = "np.sum({0},1)".format(self.evals[ctx.argumentList().expression(0)])
            elif funName=="Range":
                rFrom = self.evals[ctx.argumentList().expression(0)]
                rTo = self.evals[ctx.argumentList().expression(1)]
                rStep = self.evals[ctx.argumentList().expression(2)]
                self.evals[ctx] = "np.arange({0},{1},{2})".format(rFrom, str(float(rTo)+float(rStep)), rStep)
            elif funName=="VCat":
                self.evals[ctx] = "ml.vstack(({0}))".format(",".join(["ml.matrix({0})".format(self.evals[e]) for e in ctx.argumentList().expression()]))
            elif funName=="Zeros":
                if len(ctx.argumentList().expression())>1:
                    self.evals[ctx] = "ml.zeros(("+self.evals[ctx.argumentList()]+"))"
                else:
                    self.evals[ctx] = "np.zeros("+self.evals[ctx.argumentList()]+")"
            elif funName=="Flatten":
                self.evals[ctx] = self.evals[ctx.argumentList()]+".A.flatten()"
            elif funName=="LoadTrace":
                self.evals[ctx] = "np.loadtxt('{0}')".format(scriptDir+"/data/"+self.evals[ctx.argumentList()][1:-1])
            elif funName=="For":
                text = "for {0} in range({1},{2},{3}):".format(self.evals[ctx.argumentList().expression(0)], str(int(self.evals[ctx.argumentList().expression(1)])), str(int(self.evals[ctx.argumentList().expression(2)])+1),str(int(self.evals[ctx.argumentList().expression(3)])))
                self.forBuffer = ' '*(4*self.forDepth) + text
                self.evals[ctx] = ""               
                self.forDepth += 1
            elif funName=="EndFor":
                self.evals[ctx] = ""
                self.forDepth -= 1
                if self.forDepth==0:
                    self.evals[ctx] = self.forBuffer
                    self.forBuffer = ""
                else:
                    self.evals[ctx] = ""
            else:
                if funName in self.funMap:
                    funName = self.funMap[funName]            
                if ctx.argumentList()!=None:    
                    self.evals[ctx] = funName+"("+self.evals[ctx.argumentList()]+")"
                else:
                    self.evals[ctx] = funName+"()"
        else:
            self.evals[ctx] = "alma"

    def exitArrindexexpr (self, ctx):
        if ctx.arrix!=None:
            if ctx.arrix in self.nums:
                if ctx.arrixm!=None:
                    self.evals[ctx] = str(-self.nums[ctx.arrix]-1)
                else:
                    self.evals[ctx] = str(self.nums[ctx.arrix]-1)
            else:
                self.evals[ctx] = self.evals[ctx.arrix]+"-1"
        else:
            if ctx.arrfrom==None:
                tmpf = ""
            else:
                if ctx.arrfrom in self.nums:
                    tmpf = str(self.nums[ctx.arrfrom]-1)
                else:
                    tmpf = self.evals[ctx.arrfrom]+"-1"
            if ctx.arrende!=None:
                tmpe = self.evals[ctx.arrende]
            elif ctx.arrend==None:
                tmpe = ""
            elif ctx.arrendm!=None:
                if int(ctx.arrend.text)==0:
                    tmpe = ""
                else:
                    tmpe = "-" + str(int(ctx.arrend.text))
            else:
                tmpe = str(int(ctx.arrend.text))
                
            if ctx.arrstep==None:
                self.evals[ctx] = tmpf + ':' + tmpe
            else:
                self.evals[ctx] = tmpf + ':' + tmpe + ':' + self.evals[ctx.arrstep]

    def exitPrimary (self, ctx):
        pText = ctx.getChild(0).getText()
        if pText=="Infinity":
            pText = "np.inf"
        self.evals[ctx] = pText

# ==========================================================================================
class MATLABWriter(TestWriterBase):

    funMap = {"Norm":"norm", "Abs":"abs", "Max":"max", "Min":"min", "Total":"sum",\
              "Diff":"diff", "ToArray":"np.array", "Eigenvalues":"eig", \
              "BuTools`CheckPrecision":"BuToolsCheckPrecision", "Mean":"mean", "Inv":"inv",\
              "ExpMat":"expm","Sum":"sum","Ceil":"ceil","Linspace":"linspace",\
              "Eye":"eye","BuTools`Verbose":"BuToolsVerbose",\
              "Kron":"kron","Diag":"diag","Pinv":"pinv","MatrixMin":"np.min",\
              "MatrixMax":"np.max","Length":"length"}
    before = ""
    after = ""
    indent = 0
    printCmd = "disp({0})"
    strPrintCmd = "disp('{0}')"
    stringDelim = "'"
    vectorPre = "["
    vectorPost = "]"
    matrixPre = "["
    matrixPost = "]"
    matrixColBreak = "; "
    listPre = "["
    listPost = "]"
    cellPre = "{"
    cellPost = "}"
    rangeCmd = "({0}:{2}:{1})"
    assertCmd = "assert({0}, '{1}');"
    plotCmd =   "plot({0})"
    listExtract = "{:}"
    lineDelim = ";"
    unPackOutPre = "["
    unPackOutPost = "]"

    def printQuotedMsg(self, what):
        what = re.compile("'").sub("''", what)
        what = re.compile('"').sub("'", what)
        self.printLine("disp('{0}');".format(what))

    def printMsg(self, what):
        what = re.compile("'").sub("''", what)
        what = re.compile('"').sub("'", what)
        text = "disp({0});".format(what)
        if forDocEx:
            self.printQuotedMsg(">>> " + text)
        self.printLine(text)

    def enterTestfile (self, ctx):
        self.packName = str(ctx.ID())
        self.indent = 0
        btDir = os.path.realpath(os.path.dirname(os.path.realpath(__file__))+'/../Matlab')
        print("clear all")
        print("run('{0}/BuToolsInit.m')".format(btDir))
        print(self.before.format(self.packName))
        self.printQuotedMsg("---BuTools: {0} package test file---".format(self.packName))
        self.printQuotedMsg("Enable the verbose messages with the BuToolsVerbose flag")
        print(' '*self.indent + "global BuToolsVerbose;")
        print(' '*self.indent + "BuToolsVerbose = true;")
        self.printQuotedMsg("Enable input parameter checking with the BuToolsCheckInput flag")
        print(' '*self.indent + "global BuToolsCheckInput;")
        print(' '*self.indent + "BuToolsCheckInput = true;")
        print(' '*self.indent + "global BuToolsCheckPrecision;")
        print(' '*self.indent + "format short g;")
        if forDocEx:
            docExDir = os.path.realpath(os.path.dirname(os.path.realpath(__file__))+'/docex')+"/"
            print(' '*self.indent + "format compact")
            print(' '*self.indent + "delete('{0}_matlab.docex');".format(docExDir+self.packName))
            print(' '*self.indent + "diary('{0}_matlab.docex');".format(docExDir+self.packName))

    def exitTestfile (self, ctx):
        print(self.after.format(self.packName))

    def exitTprint (self, ctx):
        announce = str(ctx.ANNOUNCE())=="!"
        for e in ctx.expression():
            s = self.evals[e]
            if len(s)>0 and s[0]!="'" and announce and not forDocEx:
                self.printQuotedMsg(s+":")
            if s[0]!="'" or not forDocEx:
                if s[0]=="'":
                    s = s[1:-1]
                    self.printQuotedMsg(s)
                else:
                    self.printMsg(s)

    def exitExpression (self, ctx):
        if ctx.primary()!=None:
            self.evals[ctx] = self.evals[ctx.primary()]            
        elif ctx.biop!=None:
            if ctx.biop.text==".":
                oper = "*"
            else:
                oper = ctx.biop.text
            self.evals[ctx] = self.evals[ctx.expression(0)]+oper+ self.evals[ctx.expression(1)]
        elif ctx.unop!=None:
            self.evals[ctx] = ctx.unop.text + self.evals[ctx.expression(0)]
        elif ctx.par1!=None:
            self.evals[ctx] = "(" + self.evals[ctx.expression(0)]+")"
        elif ctx.arrix!=None:
            # check if we are indexing the return value of a function
            if ctx.expression(0).funpar!=None:
                if len(ctx.arrindexexpr())==1 and self.evals[ctx.arrindexexpr(0)]=="1":
                    self.evals[ctx] = self.evals[ctx.expression(0)]
                else:
                    print("function indexing!!!",file=sys.stderr)
            else:
                self.evals[ctx] = self.evals[ctx.expression(0)] + "(" + ", ".join([self.evals[e] for e in ctx.arrindexexpr()]) + ")"
        elif ctx.lstix!=None:
            self.evals[ctx] = self.evals[ctx.expression(0)] + "{" + ", ".join([self.evals[e] for e in ctx.arrindexexpr()]) + "}"
        elif ctx.tlist()!=None:
            self.evals[ctx] = self.evals[ctx.tlist()]
        elif ctx.funpar!=None:
            # function call!!!
            funName = ctx.ID().getText()
            if funName=="AllNonNegative":
                self.evals[ctx] = "all({0}>=0)".format(self.evals[ctx.argumentList()])
            elif funName=="AllNegative":
                self.evals[ctx] = "all({0}<0)".format(self.evals[ctx.argumentList()])
            elif funName=="AllGreaterThanm1":
                self.evals[ctx] = "all({0}>=-1)".format(self.evals[ctx.argumentList()])
            elif funName=="AllLessThan1":
                self.evals[ctx] = "all({0}<=1)".format(self.evals[ctx.argumentList()])
            elif funName=="AllPositiveMatrix":
                self.evals[ctx] = "all(all({0}>0))".format(self.evals[ctx.argumentList()])
            elif funName=="AllNonNegativeMatrix":
                self.evals[ctx] = "all(all({0}>=0))".format(self.evals[ctx.argumentList()])
            elif funName=="AllNotGreaterThan1Matrix":
                self.evals[ctx] = "all(all({0}<=1))".format(self.evals[ctx.argumentList()])
            elif funName=="Dot":
                self.evals[ctx] = "dot({0})".format(self.evals[ctx.argumentList()])
            elif funName=="Dim1":
                self.evals[ctx] = "size({0},1)".format(self.evals[ctx.argumentList().expression(0)])
            elif funName=="RowSum":
                self.evals[ctx] = "sum({0},2)".format(self.evals[ctx.argumentList().expression(0)])
            elif funName=="Range":
                rFrom = self.evals[ctx.argumentList().expression(0)]
                rTo = self.evals[ctx.argumentList().expression(1)]
                rStep = self.evals[ctx.argumentList().expression(2)]
                self.evals[ctx] = "({0}:{1}:{2})".format(rFrom, rStep, rTo)
            elif funName=="VCat":
                self.evals[ctx] = "[{0}]".format(";".join([self.evals[e] for e in ctx.argumentList().expression()]))
            elif funName=="Zeros":
                if len(ctx.argumentList().expression())>1:
                    self.evals[ctx] = "zeros("+self.evals[ctx.argumentList()]+")"
                else:
                    self.evals[ctx] = "zeros(1,"+self.evals[ctx.argumentList()]+")"
            elif funName=="Flatten":
                self.evals[ctx] = "reshape({0},1,numel({0}))".format(self.evals[ctx.argumentList()])
            elif funName=="ToArray":
                self.evals[ctx] = self.evals[ctx.argumentList()]
            elif funName=="MatrixMin":
                self.evals[ctx] = "min(min({0}))".format(self.evals[ctx.argumentList()])
            elif funName=="MatrixMax":
                self.evals[ctx] = "max(max({0}))".format(self.evals[ctx.argumentList()])
            elif funName=="LoadTrace":
                self.evals[ctx] = "dlmread('{0}')".format(scriptDir+"/data/"+self.evals[ctx.argumentList()][1:-1])
            elif funName=="For":
                text = "for {0}={1}:{2}:{3}".format(self.evals[ctx.argumentList().expression(0)], str(int(self.evals[ctx.argumentList().expression(1)])), str(int(self.evals[ctx.argumentList().expression(3)])),str(int(self.evals[ctx.argumentList().expression(2)])))
                self.forBuffer = ' '*(4*self.forDepth) + text
                self.evals[ctx] = ""               
                self.forDepth += 1
            elif funName=="EndFor":
                self.evals[ctx] = ""
                self.forDepth -= 1
                self.forBuffer += '\n' + ' '*(4*self.forDepth) + "end"
                if self.forDepth==0:
                    self.evals[ctx] = self.forBuffer
                    self.forBuffer = ""
                else:
                    self.evals[ctx] = ""
            else:
                if funName in self.funMap:
                    funName = self.funMap[funName]            
                if ctx.argumentList()!=None:    
                    self.evals[ctx] = funName+"("+self.evals[ctx.argumentList()]+")"
                else:
                    self.evals[ctx] = funName+"()"
        else:
            self.evals[ctx] = "alma"

    def exitArrindexexpr (self, ctx):
        if ctx.arrix!=None:
            if ctx.arrix in self.nums:
                if ctx.arrixm!=None:
                    if self.nums[ctx.arrix]==0:
                        self.evals[ctx] = "end"
                    else:
                        self.evals[ctx] = "end"+str(-self.nums[ctx.arrix])
                else:
                    self.evals[ctx] = str(self.nums[ctx.arrix])
            else:
                self.evals[ctx] = self.evals[ctx.arrix]
        else:
            if ctx.arrfrom==None:
                tmpf = "1"
            else:
                if ctx.arrfrom in self.nums:
                    tmpf = str(self.nums[ctx.arrfrom])
                else:
                    tmpf = self.evals[ctx.arrfrom]
            if ctx.arrende!=None:
                tmpe = self.evals[ctx.arrende]
            elif ctx.arrend==None:
                tmpe = "end"
            elif ctx.arrendm!=None:
                if int(ctx.arrend.text)==0:
                    tmpe = "end"
                else:
                    tmpe = "end-" + str(int(ctx.arrend.text))
            else:
                tmpe = str(int(ctx.arrend.text))
                
            if ctx.arrstep==None:
                self.evals[ctx] = tmpf + ':' + tmpe
            else:
                self.evals[ctx] = tmpf + ':' + self.evals[ctx.arrstep] + ':' + tmpe

    def exitPrimary (self, ctx):
        pText = ctx.getChild(0).getText()
        if pText=="Infinity":
            pText = "inf"
        elif pText=="True":
            pText = "true"
        elif pText=="False":
            pText = "false"
        elif pText[0]=='"':
            pText = "'"+pText[1:-1]+"'"
        self.evals[ctx] = pText


# ==========================================================================================
class MathematicaWriter(TestWriterBase):

    funMap = {"Inv":"Inverse","ExpMat":"MatrixExponential","Sum":"Total","Ceil":"Ceiling",\
              "Eye":"IdentityMatrix","Kron":"KroneckerProduct",\
              "Diag":"DiagonalMatrix","Pinv":"PseudoInverse","MatrixMin":"Min",\
              "MatrixMax":"Max","VCat":"Join","Diff":"Differences"}

    indent = 0
    printCmd = "Print[{0}]"
    strPrintCmd = 'Print["{0}"]'
    stringDelim = '"'
    vectorPre = "{"
    vectorPost = "}"
    matrixPre = "{{"
    matrixPost = "}}"
    matrixColBreak = "},{"
    listPre = "{"
    listPost = "}"
    cellPre = "{"
    cellPost = "}"
    rangeCmd = "Range[{0},{1},{2}]"
    assertCmd = 'Assert[{0}, "{1}"];'
    plotCmd =   "plot({0})"
    listExtract = ""
    lineDelim = ";"
    unPackOutPre = "{"
    unPackOutPost = "}"
    forStack = []


    def printQuotedMsg(self, what):
        what = re.compile('"').sub('\\"', what)
        self.printLine('Print["{0}"//OutputForm];'.format(what))

    def printMsg(self, what):
        text = "Print[{0}];".format(what)
        if forDocEx:
            self.printQuotedMsg(">>> " + text)
        self.printLine(text)

    def enterTestfile (self, ctx):
        self.packName = str(ctx.ID())
        self.indent = 0
        btDir = os.path.realpath(os.path.dirname(os.path.realpath(__file__))+'/../Mathematica')
        print('ClearAll["Global`*"]')
        print(' '*self.indent + 'AppendTo[$Path,"{0}"];'.format(btDir))
        print(' '*self.indent + "<<BuTools`")       

        self.printQuotedMsg("---BuTools: {0} package test file---".format(self.packName))
        self.printQuotedMsg("Enable the verbose messages with the BuToolsVerbose flag")
        print(' '*self.indent + "BuTools`Verbose = True;")
        self.printQuotedMsg("Enable input parameter checking with the BuToolsCheckInput flag")
        print(' '*self.indent + "BuTools`CheckInput = true;")
        print(' '*self.indent + "On[Assert];")
        if forDocEx:
            docExDir = os.path.realpath(os.path.dirname(os.path.realpath(__file__))+'/docex')+"/"
            print(' '*self.indent + "tmpOut = $Output;")
            print(' '*self.indent + 'stream = OpenWrite["{0}_mathematica.docex", FormatType -> InputForm, PageWidth -> Infinity];'.format(docExDir+self.packName))
            print(' '*self.indent + '$Output = {stream};')
            print(' '*self.indent + 'Unprotect[Print];')
            print(' '*self.indent + 'Print[args___] := Block[{$inMsg = True, result, str},')
            print(' '*self.indent + '   If[MatrixQ[args],')
            print(' '*self.indent + '       str = "{";')
            print(' '*self.indent + '       Do[str = StringJoin[str, ToString[args[[r]], FormatType -> InputForm]]; ')
            print(' '*self.indent + '            If[r < Length[args], str = StringJoin[str, ",\\n "]], {r, Length[args]}];')
            print(' '*self.indent + '            str = StringJoin[str, "}"];') 
            print(' '*self.indent + '            Print[str//OutputForm],')
            print(' '*self.indent + '            result = Print[args],')
            print(' '*self.indent + '            result = Print[args]')
            print(' '*self.indent + '        ]] /; ! TrueQ[$inMsg];')

    def exitTestfile (self, ctx):
        if forDocEx:
            print(' '*self.indent + '$Output = tmpOut;Close[stream];\n')

    def exitExpression (self, ctx):
        if ctx.primary()!=None:
            self.evals[ctx] = self.evals[ctx.primary()]            
        elif ctx.biop!=None:
            oper = ctx.biop.text
            if oper==".*":
                oper="*"
            elif oper==".^":
                oper="^"
            elif oper=="./":
                oper="/"
            self.evals[ctx] = self.evals[ctx.expression(0)]+oper+ self.evals[ctx.expression(1)]
        elif ctx.unop!=None:
            self.evals[ctx] = ctx.unop.text + self.evals[ctx.expression(0)]
        elif ctx.par1!=None:
            self.evals[ctx] = "(" + self.evals[ctx.expression(0)]+")"
        elif ctx.arrix!=None:
            self.evals[ctx] = self.evals[ctx.expression(0)] + "[[" + ", ".join([self.evals[e] for e in ctx.arrindexexpr()]) + "]]"
        elif ctx.lstix!=None:
            self.evals[ctx] = self.evals[ctx.expression(0)] + "[[" + ", ".join([self.evals[e] for e in ctx.arrindexexpr()]) + "]]"
        elif ctx.tlist()!=None:
            self.evals[ctx] = self.evals[ctx.tlist()]
        elif ctx.funpar!=None:
            # function call!!!
            funName = ctx.ID().getText()
            if funName=="AllNonNegative":
                self.evals[ctx] = "AllTrue[{0},#>=-10^-14 &]".format(self.evals[ctx.argumentList()])
            elif funName=="AllNegative":
                self.evals[ctx] = "AllTrue[{0},#<0 &]".format(self.evals[ctx.argumentList()])
            elif funName=="AllGreaterThanm1":
                self.evals[ctx] = "AllTrue[{0},#>=-1 &]".format(self.evals[ctx.argumentList()])
            elif funName=="AllLessThan1":
                self.evals[ctx] = "AllTrue[{0},#<=1 &]".format(self.evals[ctx.argumentList()])
            elif funName=="AllPositiveMatrix":
                self.evals[ctx] = "AllTrue[Flatten[{0}],#>0 &]".format(self.evals[ctx.argumentList()])
            elif funName=="AllNonNegativeMatrix":
                self.evals[ctx] = "AllTrue[Flatten[{0}],#>=0 &]".format(self.evals[ctx.argumentList()])
            elif funName=="AllNotGreaterThan1Matrix":
                self.evals[ctx] = "AllTrue[Flatten[{0}],#<=1 &]".format(self.evals[ctx.argumentList()])
            elif funName=="Dim1":
                self.evals[ctx] = "Dimensions[{0}][[1]]".format(self.evals[ctx.argumentList().expression(0)])
            elif funName=="RowSum":
                self.evals[ctx] = "Total[{0},{{2}}]".format(self.evals[ctx.argumentList().expression(0)])
            elif funName=="Range":
                rFrom = self.evals[ctx.argumentList().expression(0)]
                rTo = self.evals[ctx.argumentList().expression(1)]
                rStep = self.evals[ctx.argumentList().expression(2)]
                self.evals[ctx] = "Range[{0},{2},{1}]".format(rFrom, rStep, rTo)
            elif funName=="Zeros":
                if len(ctx.argumentList().expression())==2:
                    self.evals[ctx] = "Table[0,{{{0}}},{{{1}}}]".format(self.evals[ctx.argumentList().expression(0)], self.evals[ctx.argumentList().expression(1)])
                else:
                    self.evals[ctx] = "Table[0,{{{0}}}]".format(self.evals[ctx.argumentList()])
            elif funName=="Flatten":
                self.evals[ctx] = "Flatten[{0}]".format(self.evals[ctx.argumentList()])
            elif funName=="ToArray":
                self.evals[ctx] = self.evals[ctx.argumentList()]
            elif funName=="LoadTrace":
                self.evals[ctx] = 'Flatten[Import["{0}","CSV"]]'.format(scriptDir+"/data/"+self.evals[ctx.argumentList()][1:-1])
            elif funName=="Linspace":
                rFrom = self.evals[ctx.argumentList().expression(0)]
                rTo = self.evals[ctx.argumentList().expression(1)]
                numOfSteps = self.evals[ctx.argumentList().expression(2)]
                self.evals[ctx] = "Array[# &, {1}, {{{0}, {2}}}]".format(rFrom, numOfSteps, rTo)
            elif funName=="For":
                text = "Do["
                self.forStack.append(", {{{0},{1},{3},{2}}}];".format(self.evals[ctx.argumentList().expression(0)], str(int(self.evals[ctx.argumentList().expression(1)])), str(int(self.evals[ctx.argumentList().expression(3)])),str(int(self.evals[ctx.argumentList().expression(2)]))))
                self.forBuffer = ' '*(4*self.forDepth) + text
                self.evals[ctx] = ""               
                self.forDepth += 1
            elif funName=="EndFor":
                self.evals[ctx] = ""
                self.forDepth -= 1
                self.forBuffer += '\n' + ' '*(4*self.forDepth) + self.forStack[-1]
                self.forStack.pop()
                if self.forDepth==0:
                    self.evals[ctx] = self.forBuffer
                    self.forBuffer = ""
                else:
                    self.evals[ctx] = ""
            else:
                if funName in self.funMap:
                    funName = self.funMap[funName]            
                if ctx.argumentList()!=None:    
                    self.evals[ctx] = funName+"["+self.evals[ctx.argumentList()]+"]"
                else:
                    self.evals[ctx] = funName+"[]"
        else:
            self.evals[ctx] = "alma"

    def exitArrindexexpr (self, ctx):
        if ctx.arrix!=None:
            if ctx.arrix in self.nums:
                if ctx.arrixm!=None:
                    if self.nums[ctx.arrix]==0:
                        self.evals[ctx] = "-1"
                    else:
                        self.evals[ctx] = str(-self.nums[ctx.arrix]-1)
                else:
                    self.evals[ctx] = str(self.nums[ctx.arrix])
            else:
                self.evals[ctx] = self.evals[ctx.arrix]
        else:
            if ctx.arrfrom==None:
                tmpf = "1"
            else:
                if ctx.arrfrom in self.nums:
                    tmpf = str(self.nums[ctx.arrfrom])
                else:
                    tmpf = self.evals[ctx.arrfrom]
            if ctx.arrende!=None:
                tmpe = self.evals[ctx.arrende]
            elif ctx.arrend==None:
                tmpe = "-1"
            elif ctx.arrendm!=None:
                if int(ctx.arrend.text)==0:
                    tmpe = "-1"
                else:
                    tmpe = str(-int(ctx.arrend.text)-1)
            else:
                tmpe = str(int(ctx.arrend.text))
                
            if ctx.arrstep==None:
                self.evals[ctx] = tmpf + ';;' + tmpe
            else:
                self.evals[ctx] = tmpf + ';;' + tmpe + ';;' + self.evals[ctx.arrstep]

    def exitPrimary (self, ctx):
        pText = ctx.getChild(0).getText()
        self.evals[ctx] = pText

    def enterTplot (self, ctx):
        line = "ListLinePlot[{"
        for i in range(0,len(ctx.ID()),2):
            line += "Transpose[{{{0}, {1}}}]".format(ctx.ID(i),ctx.ID(i+1))
            if i+1<len(ctx.ID())-1:
                line += ","
        line += "}]"
        if forDocEx and not self.docExEnds:
            self.printQuotedMsg(">>> " + line)
        else:
            self.printLine(line)

 # ==========================================================================================
forDocEx = False
scriptDir = ""

def processFile(testFileName):

    # script directory
    global scriptDir
    scriptDir = dirname(realpath(__file__))

    # open test file
    input = FileStream(testFileName)
    baseName = splitext(basename(testFileName))[0]

    # parse test file
    lexer = bttestLexer(input)
    stream = CommonTokenStream(lexer)
    parser = bttestParser(stream)
    tree = parser.testfile()
    pythonWriter = PythonWriter()
    matlabWriter = MATLABWriter()
    mathematicaWriter = MathematicaWriter()
    walker = ParseTreeWalker()

    # generate package tests first
    global forDocEx
    forDocEx = False
    # python output
    outFile = open("{0}/testpackages/test_{1}.py".format(scriptDir, baseName),'w')
    with redirect_stdout(outFile):
        walker.walk(pythonWriter, tree)
    # matlab output        
    outFile = open("{0}/testpackages/test_{1}.m".format(scriptDir, baseName),'w')
    with redirect_stdout(outFile):
        walker.walk(matlabWriter, tree)
    # mathematica output        
    outFile = open("{0}/testpackages/test_{1}.ma".format(scriptDir, baseName),'w')
    with redirect_stdout(outFile):
        walker.walk(mathematicaWriter, tree)   

    # produce generator scripts for the examples of the documentation
    forDocEx = True
    # python output
    outFile = open("{0}/gendocex/gendocex_{1}.py".format(scriptDir, baseName),'w')
    with redirect_stdout(outFile):
        walker.walk(pythonWriter, tree)
    # matlab output        
    outFile = open("{0}/gendocex/gendocex_{1}.m".format(scriptDir, baseName),'w')
    with redirect_stdout(outFile):
        walker.walk(matlabWriter, tree)
    # mathematica output        
    outFile = open("{0}/gendocex/gendocex_{1}.ma".format(scriptDir, baseName),'w')
    with redirect_stdout(outFile):
        walker.walk(mathematicaWriter, tree)   

def main(argv):

    if len(argv)<2:
        print("Command line parameter missing! Please provide the test file name, or 'all' to process all test files in the current directory!")
    elif argv[1]=="all":
        for fName in os.listdir('.'):
            if fName.endswith('.test'):
                processFile(fName)
    else:
        processFile(argv[1])
 
if __name__ == '__main__':
    main(sys.argv)

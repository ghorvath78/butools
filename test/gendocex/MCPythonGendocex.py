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
from contextlib import redirect_stdout

def MCPythonGendocex():
	print('---BuTools: MC package test file---')
	print('Enable the verbose messages with the BuToolsVerbose flag')
	butools.verbose = True
	print('Enable input parameter checking with the BuToolsCheckInput flag')
	butools.checkInput = True
	np.set_printoptions(precision=5,linewidth=1024)
	outFile = open('MC_python.docex','w')
	with redirect_stdout(outFile):
		print('=== CRPSolve ===')
		print('>>> Q = ml.matrix([[-4.3, 3.5, 0.8],[-8.4, 6.5, 1.9],[17.3, -12.7, -4.6]])')
		Q = ml.matrix([[-4.3, 3.5, 0.8],[-8.4, 6.5, 1.9],[17.3, -12.7, -4.6]])
		print('>>> ret = CRPSolve(Q)')
		ret = CRPSolve(Q)
		print('>>> print(ret)')
		print(ret)
		print('>>> print(ret*Q)')
		print(ret*Q)
		print('=== DRPSolve ===')
		print('>>> Q = ml.matrix([[-0.9, 0.5, 1.4],[0.9, -0.9, 1],[0.3, 1.3, -0.6]])')
		Q = ml.matrix([[-0.9, 0.5, 1.4],[0.9, -0.9, 1],[0.3, 1.3, -0.6]])
		print('>>> ret = DRPSolve(Q)')
		ret = DRPSolve(Q)
		print('>>> print(ret)')
		print(ret)
		print('>>> print(ret*Q - ret)')
		print(ret*Q - ret)
		print('=== CTMCSolve ===')
		print('>>> Q = ml.matrix([[-0.9, 0.5, 0.4],[0.9, -0.9, 0],[0.3, 0.3, -0.6]])')
		Q = ml.matrix([[-0.9, 0.5, 0.4],[0.9, -0.9, 0],[0.3, 0.3, -0.6]])
		print('>>> ret = CTMCSolve(Q)')
		ret = CTMCSolve(Q)
		print('>>> print(ret)')
		print(ret)
		print('>>> print(ret*Q)')
		print(ret*Q)
		print('=== DTMCSolve ===')
		print('>>> Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.3, 0.4]])')
		Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.3, 0.4]])
		print('>>> ret = DTMCSolve(Q)')
		ret = DTMCSolve(Q)
		print('>>> print(ret)')
		print(ret)
		print('>>> print(ret*Q - ret)')
		print(ret*Q - ret)
		print('=== CheckGenerator ===')
		print('>>> Q = ml.matrix([[-0.9, 0.2, 0.4],[0, 0.9, 0.9],[0, 0.6, -0.6]])')
		Q = ml.matrix([[-0.9, 0.2, 0.4],[0, 0.9, 0.9],[0, 0.6, -0.6]])
		print('>>> flag = CheckGenerator(Q,True)')
		flag = CheckGenerator(Q,True)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[-0.9, 0.5, 0.4],[0.9, -0.9, 0],[0.3, 0.3, -0.6]])')
		Q = ml.matrix([[-0.9, 0.5, 0.4],[0.9, -0.9, 0],[0.3, 0.3, -0.6]])
		print('>>> flag = CheckGenerator(Q,True)')
		flag = CheckGenerator(Q,True)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[-0.9, 0.2, 0.4],[0.9, -0.9, 0],[0.3, 0.3, -0.6]])')
		Q = ml.matrix([[-0.9, 0.2, 0.4],[0.9, -0.9, 0],[0.3, 0.3, -0.6]])
		print('>>> flag = CheckGenerator(Q,True)')
		flag = CheckGenerator(Q,True)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[-0.9, 0.5, 0.4],[0.9, -1.1, 0],[0.3, 0.3, -0.6]])')
		Q = ml.matrix([[-0.9, 0.5, 0.4],[0.9, -1.1, 0],[0.3, 0.3, -0.6]])
		print('>>> flag = CheckGenerator(Q)')
		flag = CheckGenerator(Q)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[-0.9, 0.5, 0.4],[0.9, -0.9, 0],[0.3, 0.3, -0.6]])')
		Q = ml.matrix([[-0.9, 0.5, 0.4],[0.9, -0.9, 0],[0.3, 0.3, -0.6]])
		print('>>> flag = CheckGenerator(Q)')
		flag = CheckGenerator(Q)
		print('>>> print(flag)')
		print(flag)
		print('=== CheckProbMatrix ===')
		print('>>> Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, -0.1, 0.4]])')
		Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, -0.1, 0.4]])
		print('>>> flag = CheckProbMatrix(Q)')
		flag = CheckProbMatrix(Q)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.1, 0.4]])')
		Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.1, 0.4]])
		print('>>> flag = CheckProbMatrix(Q)')
		flag = CheckProbMatrix(Q)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.3, 0.4]])')
		Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.3, 0.4]])
		print('>>> flag = CheckProbMatrix(Q)')
		flag = CheckProbMatrix(Q)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.3, 0.4]])')
		Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.3, 0.4]])
		print('>>> flag = CheckProbMatrix(Q,True)')
		flag = CheckProbMatrix(Q,True)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.1, 0.4]])')
		Q = ml.matrix([[0.1, 0.5, 0.4],[0.9, 0.1, 0],[0.3, 0.1, 0.4]])
		print('>>> flag = CheckProbMatrix(Q,True)')
		flag = CheckProbMatrix(Q,True)
		print('>>> print(flag)')
		print(flag)
		print('=== CheckProbVector ===')
		print('>>> Q = ml.matrix([[1.1, -0.1]])')
		Q = ml.matrix([[1.1, -0.1]])
		print('>>> flag = CheckProbVector(Q)')
		flag = CheckProbVector(Q)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[1.1, 0.1]])')
		Q = ml.matrix([[1.1, 0.1]])
		print('>>> flag = CheckProbVector(Q)')
		flag = CheckProbVector(Q)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[1, 0]])')
		Q = ml.matrix([[1, 0]])
		print('>>> flag = CheckProbVector(Q)')
		flag = CheckProbVector(Q)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[0.9, -0.1]])')
		Q = ml.matrix([[0.9, -0.1]])
		print('>>> flag = CheckProbVector(Q,True)')
		flag = CheckProbVector(Q,True)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[0.9, 0.1]])')
		Q = ml.matrix([[0.9, 0.1]])
		print('>>> flag = CheckProbVector(Q,True)')
		flag = CheckProbVector(Q,True)
		print('>>> print(flag)')
		print(flag)
		print('>>> Q = ml.matrix([[0.8, 0.1]])')
		Q = ml.matrix([[0.8, 0.1]])
		print('>>> flag = CheckProbVector(Q,True)')
		flag = CheckProbVector(Q,True)
		print('>>> print(flag)')
		print(flag)
	outFile.close()

if __name__ == "__main__":
	MCPythonGendocex()

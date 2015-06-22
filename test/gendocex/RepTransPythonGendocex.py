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

def RepTransPythonGendocex():
	print('---BuTools: RepTrans package test file---')
	print('Enable the verbose messages with the BuToolsVerbose flag')
	butools.verbose = True
	print('Enable input parameter checking with the BuToolsCheckInput flag')
	butools.checkInput = True
	np.set_printoptions(precision=5,linewidth=1024)
	outFile = open('RepTrans_python.docex','w')
	with redirect_stdout(outFile):
		print('=== SimilarityMatrix ===')
		print('>>> A1m = ml.matrix([[0.2, 0.8, 0],[1.2, -0.4, 0.1],[-0.2, 0.7, 0.5]])')
		A1m = ml.matrix([[0.2, 0.8, 0],[1.2, -0.4, 0.1],[-0.2, 0.7, 0.5]])
		print('>>> T = ml.matrix([[1, 2, -4, 6],[0, 8, -9, 7],[-3, 7, 8, -2]])')
		T = ml.matrix([[1, 2, -4, 6],[0, 8, -9, 7],[-3, 7, 8, -2]])
		print('>>> A2m = la.pinv(T)*A1m*T')
		A2m = la.pinv(T)*A1m*T
		print('>>> B = SimilarityMatrix(A1m,A2m)')
		B = SimilarityMatrix(A1m,A2m)
		print('>>> err = la.norm(A1m*B-B*A2m)')
		err = la.norm(A1m*B-B*A2m)
		print('>>> print(err)')
		print(err)
		print('=== TransformToAcyclic ===')
		print('>>> A = ml.matrix([[-0.8, 0.8, 0],[0.1, -0.3, 0.1],[0.2, 0, -0.5]])')
		A = ml.matrix([[-0.8, 0.8, 0],[0.1, -0.3, 0.1],[0.2, 0, -0.5]])
		print('>>> B = TransformToAcyclic(A)')
		B = TransformToAcyclic(A)
		print('>>> print(B)')
		print(B)
		print('>>> Cm = SimilarityMatrix(A,B)')
		Cm = SimilarityMatrix(A,B)
		print('>>> err = la.norm(A*Cm-Cm*B)')
		err = la.norm(A*Cm-Cm*B)
		print('>>> print(err)')
		print(err)
		print('=== TransformToMonocyclic ===')
		print('>>> A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])')
		A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])
		print('>>> B = TransformToMonocyclic(A)')
		B = TransformToMonocyclic(A)
		print('>>> print(B)')
		print(B)
		print('>>> Cm = SimilarityMatrix(A,B)')
		Cm = SimilarityMatrix(A,B)
		print('>>> err = la.norm(A*Cm-Cm*B)')
		err = la.norm(A*Cm-Cm*B)
		print('>>> print(err)')
		print(err)
		print('=== ExtendToMarkovian ===')
		print('>>> alpha = ml.matrix([[0.2, 0.3, 0.5]])')
		alpha = ml.matrix([[0.2, 0.3, 0.5]])
		print('>>> A = ml.matrix([[-1, 0, 0],[0, -3, 0.6],[0, -0.6, -3]])')
		A = ml.matrix([[-1, 0, 0],[0, -3, 0.6],[0, -0.6, -3]])
		print('>>> B = TransformToMonocyclic(A)')
		B = TransformToMonocyclic(A)
		print('>>> print(B)')
		print(B)
		print('>>> Cm = SimilarityMatrix(A,B)')
		Cm = SimilarityMatrix(A,B)
		print('>>> beta = alpha*Cm')
		beta = alpha*Cm
		print('>>> print(beta)')
		print(beta)
		print('>>> m,M = ExtendToMarkovian(beta,B)')
		m,M = ExtendToMarkovian(beta,B)
		print('>>> print(m)')
		print(m)
		print('>>> print(M)')
		print(M)
		print('>>> Cm = SimilarityMatrix(B,M)')
		Cm = SimilarityMatrix(B,M)
		print('>>> err = la.norm(B*Cm-Cm*M)')
		err = la.norm(B*Cm-Cm*M)
		print('>>> print(err)')
		print(err)
		print('=== SimilarityMatrixForVectors ===')
		print('>>> vecA = ml.matrix([[0.0, 0.3, -1.5, 0.0]])')
		vecA = ml.matrix([[0.0, 0.3, -1.5, 0.0]])
		print('>>> vecB = ml.matrix([[1.0, 0.2, 0.0, 1.0]])')
		vecB = ml.matrix([[1.0, 0.2, 0.0, 1.0]])
		print('>>> vecA = vecA.T')
		vecA = vecA.T
		print('>>> vecB = vecB.T')
		vecB = vecB.T
		print('>>> B = SimilarityMatrixForVectors (vecA, vecB)')
		B = SimilarityMatrixForVectors (vecA, vecB)
		print('>>> print(B)')
		print(B)
		print('>>> err = la.norm(B*vecA-vecB)')
		err = la.norm(B*vecA-vecB)
		print('>>> print(err)')
		print(err)
	outFile.close()

if __name__ == "__main__":
	RepTransPythonGendocex()

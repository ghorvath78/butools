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
from contextlib import redirect_stdout

def TestRepTransPackage():
	print('---BuTools: RepTrans package test file---')
	print('Enable the verbose messages with the BuToolsVerbose flag')
	butools.verbose = True
	print('Enable input parameter checking with the BuToolsCheckInput flag')
	butools.checkInput = True
	np.set_printoptions(precision=5,linewidth=1024)
	print("========================================")
	print('Testing BuTools function SimilarityMatrix')
	A1m = ml.matrix([[0.2, 0.8, 0],[1.2, -0.4, 0.1],[-0.2, 0.7, 0.5]])
	print('A1m = ')
	print(A1m)
	T = ml.matrix([[1, 2, -4, 6],[0, 8, -9, 7],[-3, 7, 8, -2]])
	print('T = ')
	print(T)
	print('A2m = la.pinv(T)*A1m*T:')
	A2m = la.pinv(T)*A1m*T
	print('Test:')
	print('-----')
	print('B = SimilarityMatrix(A1m,A2m):')
	B = SimilarityMatrix(A1m,A2m)
	print('err = la.norm(A1m*B-B*A2m):')
	err = la.norm(A1m*B-B*A2m)
	print('err = ')
	print(err)
	assert err<10**-7 , 'The resulting matrix T does not satisfy A1m*T = T*A2m!'
	print("========================================")
	print('Testing BuTools function TransformToAcyclic')
	print('Input:')
	print('------')
	A = ml.matrix([[-0.8, 0.8, 0],[0.1, -0.3, 0.1],[0.2, 0, -0.5]])
	print('A = ')
	print(A)
	print('Test:')
	print('-----')
	print('B = TransformToAcyclic(A):')
	B = TransformToAcyclic(A)
	print('B = ')
	print(B)
	print('Cm = SimilarityMatrix(A,B):')
	Cm = SimilarityMatrix(A,B)
	print('err = la.norm(A*Cm-Cm*B):')
	err = la.norm(A*Cm-Cm*B)
	print('err = ')
	print(err)
	assert err<10**-7 , 'The original and the transformed matrix are not similar!'
	print("========================================")
	print('Testing BuTools function TransformToMonocyclic')
	A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])
	print('A = ')
	print(A)
	print('Test:')
	print('-----')
	print('B = TransformToMonocyclic(A):')
	B = TransformToMonocyclic(A)
	print('B = ')
	print(B)
	print('Cm = SimilarityMatrix(A,B):')
	Cm = SimilarityMatrix(A,B)
	print('err = la.norm(A*Cm-Cm*B):')
	err = la.norm(A*Cm-Cm*B)
	print('err = ')
	print(err)
	assert err<10**-7 , 'The original and the transformed matrix are not similar!'
	print("========================================")
	print('Testing BuTools function ExtendToMarkovian')
	print('Input:')
	print('------')
	alpha = ml.matrix([[0.2, 0.3, 0.5]])
	A = ml.matrix([[-1, 0, 0],[0, -3, 0.6],[0, -0.6, -3]])
	B = TransformToMonocyclic(A)
	print('B = ')
	print(B)
	Cm = SimilarityMatrix(A,B)
	beta = alpha*Cm
	print('beta = ')
	print(beta)
	print('Test:')
	print('-----')
	print('m,M = ExtendToMarkovian(beta,B):')
	m,M = ExtendToMarkovian(beta,B)
	print('m = ')
	print(m)
	print('M = ')
	print(M)
	print('Cm = SimilarityMatrix(B,M):')
	Cm = SimilarityMatrix(B,M)
	print('err = la.norm(B*Cm-Cm*M):')
	err = la.norm(B*Cm-Cm*M)
	print('err = ')
	print(err)
	assert err<10**-7 , 'The original and the transformed matrix are not similar!'
	assert np.min(m)>-10**-14 , 'The initial vector is still not Markovian!'
	print("========================================")
	print('Testing BuTools function SimilarityMatrixForVectors')
	print('Input:')
	print('------')
	vecA = ml.matrix([[0.0, 0.3, -1.5, 0.0]])
	print('vecA = ')
	print(vecA)
	vecB = ml.matrix([[1.0, 0.2, 0.0, 1.0]])
	print('vecB = ')
	print(vecB)
	print('Test:')
	print('-----')
	print('vecA = vecA.T:')
	vecA = vecA.T
	print('vecB = vecB.T:')
	vecB = vecB.T
	print('B = SimilarityMatrixForVectors (vecA, vecB):')
	B = SimilarityMatrixForVectors (vecA, vecB)
	print('B = ')
	print(B)
	print('err = la.norm(B*vecA-vecB):')
	err = la.norm(B*vecA-vecB)
	print('err = ')
	print(err)
	assert la.norm(B*vecA-vecB)<10**-14 , 'The resulting matrix T does not satisfy T*vecA = vecB!'

if __name__ == "__main__":
	TestRepTransPackage()

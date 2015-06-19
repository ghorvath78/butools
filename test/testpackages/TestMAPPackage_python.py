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

def TestMAPPackage():
	print('---BuTools: MAP package test file---')
	print('Enable the verbose messages with the BuToolsVerbose flag')
	butools.verbose = True
	print('Enable input parameter checking with the BuToolsCheckInput flag')
	butools.checkInput = True
	np.set_printoptions(precision=5,linewidth=1024)
	print("========================================")
	print('Testing BuTools function MarginalDistributionFromMAP')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-0.17, 0, 0, 0.07],[0.01, -0.78, 0.03, 0.08],[0.22, 0.17, -1.1, 0.02],[0.04, 0.12, 0, -0.42]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[0, 0.06, 0, 0.04],[0.04, 0.19, 0.21, 0.22],[0.22, 0.13, 0.15, 0.19],[0.05, 0, 0.17, 0.04]])
	print('D1 = ')
	print(D1)
	print('Test:')
	print('-----')
	print('a,A = MarginalDistributionFromMAP(D0,D1):')
	a,A = MarginalDistributionFromMAP(D0,D1)
	print('a = ')
	print(a)
	print('A = ')
	print(A)
	assert Length(a)==D0.shape[0]  and  CheckPHRepresentation(a,A) , 'MarginalDistributionFromMAP returned a wrong PH representation!'
	print("========================================")
	print('Testing BuTools function MarginalMomentsFromMAP')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-0.17, 0, 0, 0.07],[0.01, -0.78, 0.03, 0.08],[0.22, 0.17, -1.1, 0.02],[0.04, 0.12, 0, -0.42]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[0, 0.06, 0, 0.04],[0.04, 0.19, 0.21, 0.22],[0.22, 0.13, 0.15, 0.19],[0.05, 0, 0.17, 0.04]])
	print('D1 = ')
	print(D1)
	print('Test:')
	print('-----')
	print('moms = MarginalMomentsFromMAP(D0,D1):')
	moms = MarginalMomentsFromMAP(D0,D1)
	print('moms = ')
	print(moms)
	assert Length(moms)==2*D0.shape[0]-1  and  CheckMoments(moms) , 'MarginalMomentsFromMAP returned wrong moments!'
	print("========================================")
	print('Testing BuTools function MarginalDistributionFromRAP')
	print('Input:')
	print('------')
	H0 = ml.matrix([[-2, 0, 0],[0, -3, 1],[0, -1, -2]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[1.8, 0.2, 0],[0.2, 1.8, 0],[0.2, 1.8, 1]])
	print('H1 = ')
	print(H1)
	print('Test:')
	print('-----')
	print('a,A = MarginalDistributionFromRAP(H0,H1):')
	a,A = MarginalDistributionFromRAP(H0,H1)
	print('a = ')
	print(a)
	print('A = ')
	print(A)
	assert Length(a)==H0.shape[0]  and  CheckMERepresentation(a,A) , 'MarginalDistributionFromRAP returned a wrong ME representation!'
	print("========================================")
	print('Testing BuTools function MarginalMomentsFromRAP')
	print('Input:')
	print('------')
	H0 = ml.matrix([[-2., 0, 0],[0, -3., 1.],[0, -1., -2.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[1.8, 0.2, 0],[0.2, 1.8, 0],[0.2, 1.8, 1.]])
	print('H1 = ')
	print(H1)
	print('Test:')
	print('-----')
	print('moms = MarginalMomentsFromRAP(H0,H1):')
	moms = MarginalMomentsFromRAP(H0,H1)
	print('moms = ')
	print(moms)
	assert Length(moms)==2*H0.shape[0]-1  and  CheckMoments(moms) , 'MarginalMomentsFromRAP returned wrong moments!'
	print("========================================")
	print('Testing BuTools function MarginalDistributionFromMMAP')
	D0 = ml.matrix([[-1.78, 0.29],[0.07, -0.92]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[0.15, 0.49],[0.23, 0.36]])
	print('D1 = ')
	print(D1)
	D2 = ml.matrix([[0.11, 0.2],[0.01, 0]])
	print('D2 = ')
	print(D2)
	D3 = ml.matrix([[0.14, 0.4],[0.11, 0.14]])
	print('D3 = ')
	print(D3)
	print('Test:')
	print('-----')
	print('a,A = MarginalDistributionFromMMAP([D0,D1,D2,D3]):')
	a,A = MarginalDistributionFromMMAP([D0,D1,D2,D3])
	print('a = ')
	print(a)
	print('A = ')
	print(A)
	assert Length(a)==D0.shape[0]  and  CheckPHRepresentation(a,A) , 'MarginalDistributionFromMMAP returned a wrong PH representation!'
	print("========================================")
	print('Testing BuTools function MarginalMomentsFromMMAP')
	D0 = ml.matrix([[-1.78, 0.29],[0.07, -0.92]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[0.15, 0.49],[0.23, 0.36]])
	print('D1 = ')
	print(D1)
	D2 = ml.matrix([[0.11, 0.2],[0.01, 0]])
	print('D2 = ')
	print(D2)
	D3 = ml.matrix([[0.14, 0.4],[0.11, 0.14]])
	print('D3 = ')
	print(D3)
	print('Test:')
	print('-----')
	print('moms = MarginalMomentsFromMMAP([D0,D1,D2,D3]):')
	moms = MarginalMomentsFromMMAP([D0,D1,D2,D3])
	print('moms = ')
	print(moms)
	assert Length(moms)==2*D0.shape[0]-1  and  CheckMoments(moms) , 'MarginalMomentsFromMMAP returned wrong moments!'
	print("========================================")
	print('Testing BuTools function MarginalDistributionFromMRAP')
	print('Input:')
	print('------')
	x = 0.18
	H0 = ml.matrix([[-5., 0.1+x, 0.9, 1.],[1., -8., 0.9, 0.1],[0.9, 0.1, -4., 1.],[1., 2., 3., -9.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1., 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
	print('H1 = ')
	print(H1)
	H2 = ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1., 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1., 0.1, 0.8]])
	print('H2 = ')
	print(H2)
	print('Test:')
	print('-----')
	print('a,A = MarginalDistributionFromMRAP([H0,H1,H2]):')
	a,A = MarginalDistributionFromMRAP([H0,H1,H2])
	print('a = ')
	print(a)
	print('A = ')
	print(A)
	assert Length(a)==H0.shape[0]  and  CheckMERepresentation(a,A) , 'MarginalDistributionFromMRAP returned a wrong ME representation!'
	print("========================================")
	print('Testing BuTools function MarginalMomentsFromMRAP')
	print('Input:')
	print('------')
	x = 0.18
	H0 = ml.matrix([[-5., 0.1+x, 0.9, 1.],[1., -8., 0.9, 0.1],[0.9, 0.1, -4., 1.],[1., 2., 3., -9.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1., 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
	print('H1 = ')
	print(H1)
	H2 = ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1., 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1., 0.1, 0.8]])
	print('H2 = ')
	print(H2)
	print('Test:')
	print('-----')
	print('moms = MarginalMomentsFromMRAP([H0,H1,H2]):')
	moms = MarginalMomentsFromMRAP([H0,H1,H2])
	print('moms = ')
	print(moms)
	assert Length(moms)==2*H0.shape[0]-1  and  CheckMoments(moms) , 'MarginalMomentsFromMRAP returned wrong moments!'
	print("========================================")
	print('Testing BuTools function LagCorrelationsFromMAP')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-5., 0, 1., 1.],[1., -8., 1., 0],[1., 0, -4., 1.],[1., 2., 3., -9.]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[0, 1., 0, 2.],[2., 1., 3., 0],[0, 0, 1., 1.],[1., 1., 0, 1.]])
	print('D1 = ')
	print(D1)
	print('Test:')
	print('-----')
	print('corr = LagCorrelationsFromMAP(D0,D1,3):')
	corr = LagCorrelationsFromMAP(D0,D1,3)
	print('corr = ')
	print(corr)
	assert Length(corr)==3  and  np.all(np.array(corr)<1)  and  np.all(np.array(corr)>-1) , 'LagCorrelationsFromMAP returned wrong autocorrelation coefficients!'
	print("========================================")
	print('Testing BuTools function LagCorrelationsFromRAP')
	print('Input:')
	print('------')
	H0 = ml.matrix([[-2., 0, 0],[0, -3., 1.],[0, -1., -2.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[1.8, 0.2, 0],[0.2, 1.8, 0],[0.2, 1.8, 1.]])
	print('H1 = ')
	print(H1)
	print('Test:')
	print('-----')
	print('corr = LagCorrelationsFromRAP(H0,H1,3):')
	corr = LagCorrelationsFromRAP(H0,H1,3)
	print('corr = ')
	print(corr)
	assert Length(corr)==3  and  np.all(np.array(corr)<1)  and  np.all(np.array(corr)>-1) , 'LagCorrelationsFromRAP returned wrong autocorrelation coefficients!'
	print("========================================")
	print('Testing BuTools function LagkJointMomentsFromMAP')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-5., 0, 1., 1.],[1., -8., 1., 0],[1., 0, -4., 1.],[1., 2., 3., -9.]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[0, 1., 0, 2.],[2., 1., 3., 0],[0, 0, 1., 1.],[1., 1., 0, 1.]])
	print('D1 = ')
	print(D1)
	print('Test:')
	print('-----')
	print('Nm = LagkJointMomentsFromMAP(D0,D1,4,1):')
	Nm = LagkJointMomentsFromMAP(D0,D1,4,1)
	print('Nm = ')
	print(Nm)
	print('moms = MarginalMomentsFromMAP(D0,D1,4):')
	moms = MarginalMomentsFromMAP(D0,D1,4)
	print('moms = ')
	print(moms)
	cjm=np.empty(3)
	for i in range(3):
	    Nx=LagkJointMomentsFromMAP(D0,D1,1,i+1)
	    cjm[i] = (Nx[1,1]-moms[0]**2) / (moms[1]-moms[0]**2)
	print('Correlation from joint moments:')
	print(cjm)
	print('corr = LagCorrelationsFromMAP(D0,D1,3):')
	corr = LagCorrelationsFromMAP(D0,D1,3)
	print('corr = ')
	print(corr)
	mNm1, mNm2 = Nm[0,1:].A.flatten(), Nm[1:,0].A.flatten()
	assert np.all(Nm>0)  and  la.norm(np.array(moms)-mNm1)<10**-14  and  la.norm(np.array(moms)-mNm2)<10**-14  and  la.norm(corr-cjm)<10**-14 , 'Joint moment matrix is invalid!'
	print("========================================")
	print('Testing BuTools function LagkJointMomentsFromRAP')
	print('Input:')
	print('------')
	H0 = ml.matrix([[-2., 0, 0],[0, -3., 1.],[0, -1., -2.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[1.8, 0.2, 0],[0.2, 1.8, 0],[0.2, 1.8, 1.]])
	print('H1 = ')
	print(H1)
	print('Test:')
	print('-----')
	print('Nm = LagkJointMomentsFromRAP(H0,H1,4,1):')
	Nm = LagkJointMomentsFromRAP(H0,H1,4,1)
	print(Length(Nm))
	print('moms = MarginalMomentsFromRAP(H0,H1,4):')
	moms = MarginalMomentsFromRAP(H0,H1,4)
	print('moms = ')
	print(moms)
	cjm=np.empty(3)
	for i in range(3):
	    Nx=LagkJointMomentsFromRAP(H0,H1,1,i+1)
	    cjm[i] = (Nx[1,1]-moms[0]**2) / (moms[1]-moms[0]**2)
	print('Correlation from joint moments:')
	print(cjm)
	print('corr = LagCorrelationsFromRAP(H0,H1,3):')
	corr = LagCorrelationsFromRAP(H0,H1,3)
	print('corr = ')
	print(corr)
	mNm1, mNm2 = Nm[0,1:].A.flatten(), Nm[1:,0].A.flatten()
	assert np.all(Nm>0)  and  la.norm(np.array(moms)-mNm1)<10**-14  and  la.norm(np.array(moms)-mNm2)<10**-14  and  la.norm(corr-cjm)<10**-14 , 'Joint moment matrix is invalid!'
	print("========================================")
	print('Testing BuTools function LagkJointMomentsFromMMAP')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-1.78, 0.29],[0.07, -0.92]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[0.15, 0.49],[0.23, 0.36]])
	print('D1 = ')
	print(D1)
	D2 = ml.matrix([[0.11, 0.2],[0.01, 0]])
	print('D2 = ')
	print(D2)
	D3 = ml.matrix([[0.14, 0.4],[0.11, 0.14]])
	print('D3 = ')
	print(D3)
	print('Test:')
	print('-----')
	print('Nm = LagkJointMomentsFromMMAP([D0,D1,D2,D3],3,1):')
	Nm = LagkJointMomentsFromMMAP([D0,D1,D2,D3],3,1)
	print(Nm[0])
	print(Nm[1])
	print(Nm[2])
	assert Length(Nm)==3  and  la.norm(Nm[0]+Nm[1]+Nm[2]-LagkJointMomentsFromMAP(D0,D1+D2+D3,3,1))<10**-14 , 'Joint moment matrix is invalid!'
	print("========================================")
	print('Testing BuTools function LagkJointMomentsFromMRAP')
	print('Input:')
	print('------')
	x = 0.18
	H0 = ml.matrix([[-5., 0.1+x, 0.9, 1.],[1., -8., 0.9, 0.1],[0.9, 0.1, -4., 1.],[1., 2., 3., -9.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1., 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
	print('H1 = ')
	print(H1)
	H2 = ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1., 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1., 0.1, 0.8]])
	print('H2 = ')
	print(H2)
	print('Test:')
	print('-----')
	print('Nm = LagkJointMomentsFromMRAP([H0,H1,H2],3,2):')
	Nm = LagkJointMomentsFromMRAP([H0,H1,H2],3,2)
	print(Nm[0])
	print(Nm[1])
	assert Length(Nm)==2  and  la.norm(Nm[0]+Nm[1]-LagkJointMomentsFromRAP(H0,H1+H2,3,2))<10**-14 , 'Joint moment matrix is invalid!'
	print("========================================")
	print('Testing BuTools function RandomMAP')
	print('Test:')
	print('-----')
	print('D0,D1 = RandomMAP(4,1.62,10):')
	D0,D1 = RandomMAP(4,1.62,10)
	print('D0 = ')
	print(D0)
	print('D1 = ')
	print(D1)
	print('m = MarginalMomentsFromMAP(D0,D1,1)[0]:')
	m = MarginalMomentsFromMAP(D0,D1,1)[0]
	print('m = ')
	print(m)
	assert CheckMAPRepresentation(D0,D1) , 'RandomMAP failed to return a valid MAP representation!'
	assert np.abs(m-1.62)<10**-14 , 'RandomMAP failed to match the given mean value!'
	print("========================================")
	print('Testing BuTools function RandomMMAP')
	print('Test:')
	print('-----')
	print('D = RandomMMAP(4,3,1.62,10):')
	D = RandomMMAP(4,3,1.62,10)
	print(D[0])
	print(D[1])
	print(D[2])
	print(D[3])
	print('m = MarginalMomentsFromMMAP(D,1)[0]:')
	m = MarginalMomentsFromMMAP(D,1)[0]
	print('m = ')
	print(m)
	assert CheckMMAPRepresentation(D) , 'RandomMMAP failed to return a valid MMAP representation!'
	assert np.abs(m-1.62)<10**-14 , 'RandomMMAP failed to match the given mean value!'
	print("========================================")
	print('Testing BuTools function CheckMAPRepresentation')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-1., 0, 1.],[0, -2., 0],[1., 0, -3.]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[-1., 0, 1., 0],[0, -2., 0, 1.],[1., 0, -3., 0],[1., 2., 2., 1.]])
	print('D1 = ')
	print(D1)
	print('Test:')
	print('-----')
	print('flag = CheckMAPRepresentation(D0,D1):')
	flag = CheckMAPRepresentation(D0,D1)
	print('flag = ')
	print(flag)
	assert flag==False , 'CheckMAPRepresentation failed to detect the incompatible shapes of D0 and D1!'
	print('Input:')
	print('------')
	D0 = ml.matrix([[-1., 0, 1.],[0, -2., 0],[1., 0, -3.]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[1., 0, 1.],[0, 2., 0],[1., 0, 3.]])
	print('D1 = ')
	print(D1)
	print('Test:')
	print('-----')
	print('flag = CheckMAPRepresentation(D0,D1):')
	flag = CheckMAPRepresentation(D0,D1)
	print('flag = ')
	print(flag)
	assert flag==False , 'CheckMAPRepresentation failed to detect invalid rowsums!'
	print('Input:')
	print('------')
	D0 = ml.matrix([[-3., 0, 1.],[0, -2., 0],[1., 0, -5.]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[1., 0, 1.],[0, 2., 0],[1., 0, 3.]])
	print('D1 = ')
	print(D1)
	print('Test:')
	print('-----')
	print('flag = CheckMAPRepresentation(D0,D1):')
	flag = CheckMAPRepresentation(D0,D1)
	print('flag = ')
	print(flag)
	assert flag==True , 'CheckMAPRepresentation failed to recognize a valid MAP representation!'
	print("========================================")
	print('Testing BuTools function CheckRAPRepresentation')
	print('Input:')
	print('------')
	H0 = ml.matrix([[-1., 0, 1.],[0, -2., 0],[1., 0, -3.],[1., 2., 2.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[-1., 0, 1.],[0, -2., 0],[1., 0, -3.],[1., 2., 2.]])
	print('H1 = ')
	print(H1)
	print('Test:')
	print('-----')
	print('flag = CheckRAPRepresentation(H0,H1):')
	flag = CheckRAPRepresentation(H0,H1)
	print('flag = ')
	print(flag)
	assert flag==False , 'CheckRAPRepresentation failed to detect the incompatible shapes of D0 and D1!'
	print('Input:')
	print('------')
	H0 = ml.matrix([[-1., 0, 2.],[0, 2., 0],[1., 0, -3.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[-1., 0, 1.],[0, -2., 0],[1., 0, -3.]])
	print('H1 = ')
	print(H1)
	print('Test:')
	print('-----')
	print('flag = CheckRAPRepresentation(H0,H1):')
	flag = CheckRAPRepresentation(H0,H1)
	print('flag = ')
	print(flag)
	assert flag==False , 'CheckRAPRepresentation failed to detect invalid rowsums!'
	print('Input:')
	print('------')
	H0 = ml.matrix([[-1., 0, 0],[0, -2., 2.],[0, 3., -3.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[0, 0, 1.],[0, -1., 1.],[1., 0, -1.]])
	print('H1 = ')
	print(H1)
	print('Test:')
	print('-----')
	print('flag = CheckRAPRepresentation(H0,H1):')
	flag = CheckRAPRepresentation(H0,H1)
	print('flag = ')
	print(flag)
	assert flag==False , 'CheckRAPRepresentation failed to detect invalid eigenvalues!'
	print('Input:')
	print('------')
	H0 = ml.matrix([[-2., 0, 0],[0, -1., 1.],[0, -1., -1.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[1., 0, 1.],[0, 1., -1.],[1., 0, 1.]])
	print('H1 = ')
	print(H1)
	print('Test:')
	print('-----')
	print('flag = CheckRAPRepresentation(H0,H1):')
	flag = CheckRAPRepresentation(H0,H1)
	print('flag = ')
	print(flag)
	assert flag==False , 'CheckRAPRepresentation failed to detect invalid eigenvalues!'
	print('Input:')
	print('------')
	H0 = ml.matrix([[-1., 0, 0],[0, -2., 1.],[0, -1., -2.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[1., 0, 0],[0, 1., 0],[1., 1., 1.]])
	print('H1 = ')
	print(H1)
	print('Test:')
	print('-----')
	print('flag = CheckRAPRepresentation(H0,H1):')
	flag = CheckRAPRepresentation(H0,H1)
	print('flag = ')
	print(flag)
	assert flag==True , 'CheckRAPRepresentation failed to recognize a valid RAP representation!'
	print("========================================")
	print('Testing BuTools function CheckMMAPRepresentation')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-1.05, 0.03, 0.07],[0.19, -1.63, 0.06],[0, 0.2, -1.03]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[0.16, 0.11, 0],[0.1, 0.16, 0],[0.27, 0, 0.19]])
	print('D1 = ')
	print(D1)
	D2 = ml.matrix([[0.01, 0.09, 0.13],[0.26, 0.21, 0.05],[0, 0.16, 0.07]])
	print('D2 = ')
	print(D2)
	D3 = ml.matrix([[0.19, 0.06, 0.2],[0.17, 0.16, 0.27],[0, 0, 0.14]])
	print('D3 = ')
	print(D3)
	print('Test:')
	print('-----')
	print('flag = CheckMMAPRepresentation([D0,D1,D2,D3]):')
	flag = CheckMMAPRepresentation([D0,D1,D2,D3])
	print('flag = ')
	print(flag)
	assert flag==True , 'CheckMMAPRepresentation failed to recognize a valid MMAP representation!'
	print("========================================")
	print('Testing BuTools function CheckMRAPRepresentation')
	print('Input:')
	print('------')
	x = 0.18
	H0 = ml.matrix([[-5., 0.1+x, 0.9, 1.],[1., -8., 0.9, 0.1],[0.9, 0.1, -4., 1.],[1., 2., 3., -9.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1., 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
	print('H1 = ')
	print(H1)
	H2 = ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1., 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1., 0.1, 0.8]])
	print('H2 = ')
	print(H2)
	print('Test:')
	print('-----')
	print('flag = CheckMRAPRepresentation([H0,H1,H2]):')
	flag = CheckMRAPRepresentation([H0,H1,H2])
	print('flag = ')
	print(flag)
	assert flag==True , 'CheckMRAPRepresentation failed to recognize a valid MRAP representation!'
	print("========================================")
	print('Testing BuTools function RAPFromMoments')
	print('Input:')
	print('------')
	G0 = ml.matrix([[-6.2, 2., 0.],[2., -9., 1.],[1., 0, -3.]])
	print('G0 = ')
	print(G0)
	G1 = ml.matrix([[2.2, -2., 4.],[2., 2., 2.],[1., 0, 1.]])
	print('G1 = ')
	print(G1)
	print('moms = MarginalMomentsFromRAP(G0,G1,5):')
	moms = MarginalMomentsFromRAP(G0,G1,5)
	print('moms = ')
	print(moms)
	print('Nm = LagkJointMomentsFromRAP(G0,G1,2,1):')
	Nm = LagkJointMomentsFromRAP(G0,G1,2,1)
	print('Nm = ')
	print(Nm)
	print('Test:')
	print('-----')
	print('H0,H1 = RAPFromMoments(moms,Nm):')
	H0,H1 = RAPFromMoments(moms,Nm)
	print('H0 = ')
	print(H0)
	print('H1 = ')
	print(H1)
	print('rmoms = MarginalMomentsFromRAP(H0,H1,5):')
	rmoms = MarginalMomentsFromRAP(H0,H1,5)
	print('rmoms = ')
	print(rmoms)
	print('rNm = LagkJointMomentsFromRAP(H0,H1,2,1):')
	rNm = LagkJointMomentsFromRAP(H0,H1,2,1)
	print('rNm = ')
	print(rNm)
	assert la.norm(np.array(moms)-np.array(rmoms))<10**-12  and  la.norm(Nm-rNm)<10**-12 , 'The moments and joint moments returned by RAPFromMoments are not the same as given!'
	print('Input:')
	print('------')
	G0 = ml.matrix([[-5., 0, 1., 1.],[1., -8., 1., 0],[1., 0, -4., 1.],[1., 2., 3., -9.]])
	print('G0 = ')
	print(G0)
	G1 = ml.matrix([[0, 1., 0, 2.],[2., 1., 3., 0],[0, 0, 1., 1.],[1., 1., 0, 1.]])
	print('G1 = ')
	print(G1)
	print('moms = MarginalMomentsFromRAP(G0,G1,7):')
	moms = MarginalMomentsFromRAP(G0,G1,7)
	print('moms = ')
	print(moms)
	print('Nm = LagkJointMomentsFromRAP(G0,G1,3,1):')
	Nm = LagkJointMomentsFromRAP(G0,G1,3,1)
	print('Nm = ')
	print(Nm)
	print('Test:')
	print('-----')
	print('H0,H1 = RAPFromMoments(moms,Nm):')
	H0,H1 = RAPFromMoments(moms,Nm)
	print('H0 = ')
	print(H0)
	print('H1 = ')
	print(H1)
	butools.checkPrecision = 10.**-8
	print('rmoms = MarginalMomentsFromRAP(H0,H1,7):')
	rmoms = MarginalMomentsFromRAP(H0,H1,7)
	print('rmoms = ')
	print(rmoms)
	print('rNm = LagkJointMomentsFromRAP(H0,H1,3,1):')
	rNm = LagkJointMomentsFromRAP(H0,H1,3,1)
	print('rNm = ')
	print(rNm)
	assert CheckRAPRepresentation(H0,H1,10**-8) , 'RAPFromMoments returned an invalid RAP representation!'
	assert la.norm(np.array(moms)-np.array(rmoms))<10**-8  and  la.norm(Nm-rNm)<10**-8 , 'The moments and joint moments returned by RAPFromMoments are not the same as given!'
	print("========================================")
	print('Testing BuTools function MRAPFromMoments')
	print('Input:')
	print('------')
	G0 = ml.matrix([[-1.05, 0.03, 0.07],[0.19, -1.63, 0.06],[0, 0.2, -1.03]])
	print('G0 = ')
	print(G0)
	G1 = ml.matrix([[0.16, 0.11, 0],[0.1, 0.16, 0],[0.27, 0, 0.19]])
	print('G1 = ')
	print(G1)
	G2 = ml.matrix([[0.01, 0.09, 0.13],[0.26, 0.21, 0.05],[0, 0.16, 0.07]])
	print('G2 = ')
	print(G2)
	G3 = ml.matrix([[0.19, 0.06, 0.2],[0.17, 0.16, 0.27],[0, 0, 0.14]])
	print('G3 = ')
	print(G3)
	print('G = [G0,G1,G2,G3]:')
	G = [G0,G1,G2,G3]
	print('moms = MarginalMomentsFromMRAP(G,5):')
	moms = MarginalMomentsFromMRAP(G,5)
	print('moms = ')
	print(moms)
	print('Nm = LagkJointMomentsFromMRAP(G,2,1):')
	Nm = LagkJointMomentsFromMRAP(G,2,1)
	print('Nm1, Nm2, Nm3 = Nm:')
	Nm1, Nm2, Nm3 = Nm
	print('Nm1 = ')
	print(Nm1)
	print('Nm2 = ')
	print(Nm2)
	print('Nm3 = ')
	print(Nm3)
	print('Test:')
	print('-----')
	print('H = MRAPFromMoments(moms,Nm):')
	H = MRAPFromMoments(moms,Nm)
	print('H[0]:')
	print(H[0])
	print('H[1]:')
	print(H[1])
	print('H[2]:')
	print(H[2])
	print('H[3]:')
	print(H[3])
	butools.checkPrecision = 10.**-10
	print('rmoms = MarginalMomentsFromMRAP(H,5):')
	rmoms = MarginalMomentsFromMRAP(H,5)
	print('rmoms = ')
	print(rmoms)
	print('rNm = LagkJointMomentsFromMRAP(H,2,1):')
	rNm = LagkJointMomentsFromMRAP(H,2,1)
	print('rNm1, rNm2, rNm3 = rNm:')
	rNm1, rNm2, rNm3 = rNm
	print('rNm1 = ')
	print(rNm1)
	print('rNm2 = ')
	print(rNm2)
	print('rNm3 = ')
	print(rNm3)
	assert la.norm(np.array(moms)-np.array(rmoms))<10**-9  and  la.norm(Nm1-rNm1)<10**-10  and  la.norm(Nm2-rNm2)<10**-10  and  la.norm(Nm3-rNm3)<10**-10 , 'The moments and joint moments returned by MRAPFromMoments are not the same as given!'
	print("========================================")
	print('Testing BuTools function RAPFromMomentsAndCorrelations')
	print('Input:')
	print('------')
	H0 = ml.matrix([[-6.2, 2., 0],[2., -9., 1.],[1., 0, -3.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[2.2, 0, 2.],[0, 4., 2.],[0, 1., 1.]])
	print('H1 = ')
	print(H1)
	print('mom = MarginalMomentsFromRAP(H0,H1):')
	mom = MarginalMomentsFromRAP(H0,H1)
	print('mom = ')
	print(mom)
	print('corr = LagCorrelationsFromRAP(H0,H1,3):')
	corr = LagCorrelationsFromRAP(H0,H1,3)
	print('corr = ')
	print(corr)
	print('Test:')
	print('-----')
	print('G0,G1 = RAPFromMomentsAndCorrelations(mom,corr):')
	G0,G1 = RAPFromMomentsAndCorrelations(mom,corr)
	print('G0 = ')
	print(G0)
	print('G1 = ')
	print(G1)
	print('rmom = MarginalMomentsFromRAP(G0,G1):')
	rmom = MarginalMomentsFromRAP(G0,G1)
	print('rmom = ')
	print(rmom)
	print('rcorr = LagCorrelationsFromRAP(G0,G1,3):')
	rcorr = LagCorrelationsFromRAP(G0,G1,3)
	print('rcorr = ')
	print(rcorr)
	assert CheckRAPRepresentation(G0,G1) , 'RAPFromMomentsAndCorrelations returned an invalid RAP representation!'
	assert la.norm(np.array(rmom)-np.array(mom))+la.norm(np.array(rcorr)-np.array(corr))<10**-12 , 'The result of RAPFromMomentsAndCorrelations has different moments or correlations than given!'
	print("========================================")
	print('Testing BuTools function MAP2FromMoments')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-14., 1.],[1, -25.]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[6., 7.],[3., 21.]])
	print('D1 = ')
	print(D1)
	print('moms = MarginalMomentsFromMAP(D0, D1, 3):')
	moms = MarginalMomentsFromMAP(D0, D1, 3)
	print('moms = ')
	print(moms)
	print('corr = LagCorrelationsFromMAP(D0, D1, 1)[0]:')
	corr = LagCorrelationsFromMAP(D0, D1, 1)[0]
	print('corr = ')
	print(corr)
	print('Test:')
	print('-----')
	print('D0,D1 = MAP2FromMoments(moms,corr):')
	D0,D1 = MAP2FromMoments(moms,corr)
	print('D0 = ')
	print(D0)
	print('D1 = ')
	print(D1)
	print('rmoms = MarginalMomentsFromMAP(D0, D1, 3):')
	rmoms = MarginalMomentsFromMAP(D0, D1, 3)
	print('rmoms = ')
	print(rmoms)
	print('rcorr = LagCorrelationsFromMAP(D0, D1, 1)[0]:')
	rcorr = LagCorrelationsFromMAP(D0, D1, 1)[0]
	print('rcorr = ')
	print(rcorr)
	assert CheckMAPRepresentation(D0,D1) , 'MAP2FromMoments returned an invalid MAP representation!'
	assert la.norm(np.array(moms)-np.array(rmoms))<10**-12  and  np.abs(corr-rcorr)<10**-12 , 'The moments and the correlation returned by MAP2FromMoments are not the same as given!'
	print("========================================")
	print('Testing BuTools function MAP2CorrelationBounds')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-14., 1.],[1., -25.]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[6., 7.],[3., 21.]])
	print('D1 = ')
	print(D1)
	print('moms = MarginalMomentsFromMAP(D0,D1,3):')
	moms = MarginalMomentsFromMAP(D0,D1,3)
	print('moms = ')
	print(moms)
	print('Test:')
	print('-----')
	print('lb,ub = MAP2CorrelationBounds(moms):')
	lb,ub = MAP2CorrelationBounds(moms)
	print('lb = ')
	print(lb)
	print('ub = ')
	print(ub)
	assert lb<=0  and  lb>=-1  and  ub>=0  and  ub<=1 , 'Correlation bounds given by MAP2CorrelationBounds are not correct'
	print("========================================")
	print('Testing BuTools function MAPFromFewMomentsAndCorrelations')
	print('Input:')
	print('------')
	moms = [1.1, 6.05]
	print('moms = ')
	print(moms)
	corr1 = -0.17
	print('corr1 = ')
	print(corr1)
	print('Test:')
	print('-----')
	print('D0,D1 = MAPFromFewMomentsAndCorrelations(moms, corr1):')
	D0,D1 = MAPFromFewMomentsAndCorrelations(moms, corr1)
	print('D0 = ')
	print(D0)
	print('D1 = ')
	print(D1)
	print('rmoms = MarginalMomentsFromMAP(D0,D1,2):')
	rmoms = MarginalMomentsFromMAP(D0,D1,2)
	print('rmoms = ')
	print(rmoms)
	print('rcorr1 = LagCorrelationsFromMAP(D0,D1,1)[0]:')
	rcorr1 = LagCorrelationsFromMAP(D0,D1,1)[0]
	print('rcorr1 = ')
	print(rcorr1)
	assert CheckMAPRepresentation(D0,D1) , 'MAPFromFewMomentsAndCorrelations returned with a non-Markovian representation!'
	assert la.norm(np.array(rmoms)-np.array(moms))<10**-12  and  np.abs(rcorr1-corr1)<10**-12 , 'MAPFromFewMomentsAndCorrelations failed to match the marginal moments or the lag-1 autocorrelation!'
	print('Input:')
	print('------')
	moms = [1.2, 4.32, 20.]
	print('moms = ')
	print(moms)
	corr1 = -0.4
	print('corr1 = ')
	print(corr1)
	print('Test:')
	print('-----')
	print('D0,D1 = MAPFromFewMomentsAndCorrelations(moms, corr1):')
	D0,D1 = MAPFromFewMomentsAndCorrelations(moms, corr1)
	print('D0 = ')
	print(D0)
	print('D1 = ')
	print(D1)
	butools.checkPrecision = 10.**-12
	print('rmoms = MarginalMomentsFromMAP(D0,D1,3):')
	rmoms = MarginalMomentsFromMAP(D0,D1,3)
	print('rmoms = ')
	print(rmoms)
	print('rcorr1 = LagCorrelationsFromMAP(D0,D1,1)[0]:')
	rcorr1 = LagCorrelationsFromMAP(D0,D1,1)[0]
	print('rcorr1 = ')
	print(rcorr1)
	assert CheckMAPRepresentation(D0,D1,10**-13) , 'MAPFromFewMomentsAndCorrelations returned with a non-Markovian representation!'
	assert la.norm(np.array(rmoms)-np.array(moms))<10**-12  and  np.abs(rcorr1-corr1)<10**-12 , 'MAPFromFewMomentsAndCorrelations failed to match the marginal moments or the lag-1 autocorrelation!'
	print('Input:')
	print('------')
	moms = [1.2, 4.32, 20.]
	print('moms = ')
	print(moms)
	corr1 = 0.4
	print('corr1 = ')
	print(corr1)
	print('Test:')
	print('-----')
	print('D0,D1 = MAPFromFewMomentsAndCorrelations(moms, corr1):')
	D0,D1 = MAPFromFewMomentsAndCorrelations(moms, corr1)
	print('D0 = ')
	print(D0)
	print('D1 = ')
	print(D1)
	print('rmoms = MarginalMomentsFromMAP(D0,D1,3):')
	rmoms = MarginalMomentsFromMAP(D0,D1,3)
	print('rmoms = ')
	print(rmoms)
	print('rcorr1 = LagCorrelationsFromMAP(D0,D1,1)[0]:')
	rcorr1 = LagCorrelationsFromMAP(D0,D1,1)[0]
	print('rcorr1 = ')
	print(rcorr1)
	assert CheckMAPRepresentation(D0,D1,10**-13) , 'MAPFromFewMomentsAndCorrelations returned with a non-Markovian representation!'
	assert la.norm(np.array(rmoms)-np.array(moms))<10**-12  and  np.abs(rcorr1-corr1)<10**-12 , 'MAPFromFewMomentsAndCorrelations failed to match the marginal moments or the lag-1 autocorrelation!'
	print("========================================")
	print('Testing BuTools function CanonicalFromMAP2')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-14., 1.],[1., -25.]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[6., 7.],[3., 21.]])
	print('D1 = ')
	print(D1)
	print('Test:')
	print('-----')
	print('H0,H1 = CanonicalFromMAP2(D0,D1):')
	H0,H1 = CanonicalFromMAP2(D0,D1)
	print('H0 = ')
	print(H0)
	print('H1 = ')
	print(H1)
	Cm = SimilarityMatrix(H0,D0)
	err1 = la.norm(H0*Cm-Cm*D0)
	err2 = la.norm(H1*Cm-Cm*D1)
	print('Transformation errors:')
	print(np.max(err1,err2))
	assert CheckMAPRepresentation(H0,H1) , 'The result of CanonicalFromMAP2 is not a valid MAP representation!'
	assert np.max(err1,err2)<10**-12 , 'The MAP returned by CanonicalFromMAP2 is not similar to the input!'
	print("========================================")
	print('Testing BuTools function MAPFromRAP')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-2., 2.],[2., -9.]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[-2., 2.],[3., 4.]])
	print('D1 = ')
	print(D1)
	print('Test:')
	print('-----')
	print('H0,H1 = MAPFromRAP(D0,D1):')
	H0,H1 = MAPFromRAP(D0,D1)
	print('H0 = ')
	print(H0)
	print('H1 = ')
	print(H1)
	print('err = la.norm(LagkJointMomentsFromRAP(D0,D1,3,1)-LagkJointMomentsFromRAP(H0,H1,3,1)):')
	err = la.norm(LagkJointMomentsFromRAP(D0,D1,3,1)-LagkJointMomentsFromRAP(H0,H1,3,1))
	print('err = ')
	print(err)
	assert err<10**-12 , 'The RAP returned by MAPFromRAP is not similar to the input!'
	print('Input:')
	print('------')
	D0 = ml.matrix([[-2.4, 2.],[2., -9.]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[-1.6, 2.],[3., 4.]])
	print('D1 = ')
	print(D1)
	print('Test:')
	print('-----')
	print('H0,H1 = MAPFromRAP(D0,D1):')
	H0,H1 = MAPFromRAP(D0,D1)
	print('H0 = ')
	print(H0)
	print('H1 = ')
	print(H1)
	print('err = la.norm(LagkJointMomentsFromRAP(D0,D1,3,1)-LagkJointMomentsFromRAP(H0,H1,3,1)):')
	err = la.norm(LagkJointMomentsFromRAP(D0,D1,3,1)-LagkJointMomentsFromRAP(H0,H1,3,1))
	print('err = ')
	print(err)
	assert err<10**-12 , 'The MAP returned by MAPFromRAP is not similar to the input!'
	assert CheckMAPRepresentation(H0,H1) , 'The result of MAPFromRAP is not a MAP, as it should be!'
	print("========================================")
	print('Testing BuTools function MMAPFromMRAP')
	print('Input:')
	print('------')
	x = 0.18
	H0 = ml.matrix([[-5., 0.1+x, 0.9, 1.],[1., -8., 0.9, 0.1],[0.9, 0.1, -4., 1.],[1., 2., 3., -9.]])
	print('H0 = ')
	print(H0)
	H1 = ml.matrix([[0.1-x, 0.7, 0.1, 0.1],[0.1, 1., 1.8, 0.1],[0.1, 0.1, 0.1, 0.7],[0.7, 0.1, 0.1, 0.1]])
	print('H1 = ')
	print(H1)
	H2 = ml.matrix([[0.1, 0.1, 0.1, 1.7],[1.8, 0.1, 1., 0.1],[0.1, 0.1, 0.7, 0.1],[0.1, 1., 0.1, 0.8]])
	print('H2 = ')
	print(H2)
	H = [H0,H1,H2]
	print('moms = MarginalMomentsFromMRAP(H):')
	moms = MarginalMomentsFromMRAP(H)
	print('moms = ')
	print(moms)
	print('jmom = LagkJointMomentsFromMRAP(H,3,1):')
	jmom = LagkJointMomentsFromMRAP(H,3,1)
	print('Test:')
	print('-----')
	print('G = MMAPFromMRAP(H):')
	G = MMAPFromMRAP(H)
	print(G[0])
	print(G[1])
	print(G[2])
	print('rmoms = MarginalMomentsFromMMAP(G):')
	rmoms = MarginalMomentsFromMMAP(G)
	print('rmoms = ')
	print(rmoms)
	print('rjmom = LagkJointMomentsFromMMAP(G,3,1):')
	rjmom = LagkJointMomentsFromMMAP(G,3,1)
	print('err = la.norm(rjmom[0]-jmom[0]) + la.norm(rjmom[1]-jmom[1]):')
	err = la.norm(rjmom[0]-jmom[0]) + la.norm(rjmom[1]-jmom[1])
	print('err = ')
	print(err)
	assert err<10**-12 , 'The MMAP returned by MMAPFromMRAP is not similar to the input!'
	assert CheckMMAPRepresentation(G) , 'The result of MMAPFromMRAP is not a MMAP, as it should be!'
	print("========================================")
	print('Testing BuTools function MinimalRepFromRAP')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-5., 1., 0],[3., -3., 0],[1., 1., -5.]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[0, 0, 4.],[0, 0, 0],[1., 1., 1.]])
	print('D1 = ')
	print(D1)
	print('Test:')
	print('-----')
	print('H0,H1 = MinimalRepFromRAP(D0,D1,"cont"):')
	H0,H1 = MinimalRepFromRAP(D0,D1,"cont")
	print('H0 = ')
	print(H0)
	print('H1 = ')
	print(H1)
	Cm = SimilarityMatrix(H0,D0)
	err1 = la.norm(H0*Cm-Cm*D0)
	err2 = la.norm(H1*Cm-Cm*D1)
	print('Transformation errors:')
	print(np.max(err1,err2))
	assert CheckRAPRepresentation(H0,H1) , 'MinimalRepFromRAP did not return a valid RAP representation!'
	assert H0.shape[0]==3  and  np.max(err1,err2)<10**-12 , 'MinimalRepFromRAP returned a RAP which is non-similar to the input or has an unexpected size!'
	print('H0,H1 = MinimalRepFromRAP(D0,D1,"obs"):')
	H0,H1 = MinimalRepFromRAP(D0,D1,"obs")
	print('H0 = ')
	print(H0)
	print('H1 = ')
	print(H1)
	Cm = SimilarityMatrix(H0,D0)
	err1 = la.norm(H0*Cm-Cm*D0)
	err2 = la.norm(H1*Cm-Cm*D1)
	print('Transformation errors:')
	print(np.max(err1,err2))
	assert CheckRAPRepresentation(H0,H1) , 'MinimalRepFromRAP did not return a valid RAP representation!'
	assert H0.shape[0]==2  and  np.max(err1,err2)<10**-12 , 'MinimalRepFromRAP returned a RAP which is non-similar to the input or has an unexpected size!'
	print('H0,H1 = MinimalRepFromRAP(D0,D1,"obscont"):')
	H0,H1 = MinimalRepFromRAP(D0,D1,"obscont")
	print('H0 = ')
	print(H0)
	print('H1 = ')
	print(H1)
	Cm = SimilarityMatrix(H0,D0)
	err1 = la.norm(H0*Cm-Cm*D0)
	err2 = la.norm(H1*Cm-Cm*D1)
	print('Transformation errors:')
	print(np.max(err1,err2))
	assert CheckRAPRepresentation(H0,H1) , 'MinimalRepFromRAP did not return a valid RAP representation!'
	assert H0.shape[0]==2  and  np.max(err1,err2)<10**-12 , 'MinimalRepFromRAP returned a RAP which is non-similar to the input or has an '
	print("========================================")
	print('Testing BuTools function MinimalRepFromMRAP')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-5., 1., 0],[3., -3., 0],[1., 1., -5.]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[0, 0, 0.8],[0, 0, 0],[0.2, 0.2, 0.2]])
	print('D1 = ')
	print(D1)
	D2 = ml.matrix([[0, 0, 3.2],[0, 0, 0],[0.8, 0.8, 0.8]])
	print('D2 = ')
	print(D2)
	D = [D0,D1,D2]
	print('Test:')
	print('-----')
	print('H = MinimalRepFromMRAP(D,"cont"):')
	H = MinimalRepFromMRAP(D,"cont")
	print(H[0])
	print(H[1])
	print(H[2])
	Cm = SimilarityMatrix(H[0],D[0])
	err = la.norm(H[0]*Cm-Cm*D[0]) + la.norm(H[1]*Cm-Cm*D[1]) + la.norm(H[2]*Cm-Cm*D[2])
	print('Transformation errors:')
	print(err)
	assert CheckMRAPRepresentation(H) , 'MinimalRepFromMRAP did not return a valid MRAP representation!'
	assert H[1].shape[0]==3  and  err<10**-12 , 'MinimalRepFromMRAP returned a MRAP which is non-similar to the input or has an unexpected size!'
	print('H = MinimalRepFromMRAP(D,"obs"):')
	H = MinimalRepFromMRAP(D,"obs")
	print(H[0])
	print(H[1])
	print(H[2])
	Cm = SimilarityMatrix(H[0],D[0])
	err = la.norm(H[0]*Cm-Cm*D[0]) + la.norm(H[1]*Cm-Cm*D[1]) + la.norm(H[2]*Cm-Cm*D[2])
	print('Transformation errors:')
	print(err)
	assert CheckMRAPRepresentation(H) , 'MinimalRepFromMRAP did not return a valid MRAP representation!'
	assert H[1].shape[0]==2  and  err<10**-12 , 'MinimalRepFromMRAP returned a MRAP which is non-similar to the input or has an unexpected size!'
	print('H = MinimalRepFromMRAP(D,"obscont"):')
	H = MinimalRepFromMRAP(D,"obscont")
	print(H[0])
	print(H[1])
	print(H[2])
	Cm = SimilarityMatrix(H[0],D[0])
	err = la.norm(H[0]*Cm-Cm*D[0]) + la.norm(H[1]*Cm-Cm*D[1]) + la.norm(H[2]*Cm-Cm*D[2])
	print('Transformation errors:')
	print(err)
	assert CheckMRAPRepresentation(H) , 'MinimalRepFromMRAP did not return a valid MRAP representation!'
	assert H[1].shape[0]==2  and  err<10**-12 , 'MinimalRepFromMRAP returned a MRAP which is non-similar to the input or has an unexpected size!'
	print("========================================")
	print('Testing BuTools function SamplesFromMAP')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-0.17, 0, 0, 0.07],[0.01, -0.78, 0.03, 0.08],[0.22, 0.17, -1.1, 0.02],[0.04, 0.12, 0, -0.42]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[0, 0.06, 0, 0.04],[0.04, 0.19, 0.21, 0.22],[0.22, 0.13, 0.15, 0.19],[0.05, 0, 0.17, 0.04]])
	print('D1 = ')
	print(D1)
	print('Test:')
	print('-----')
	print('x = SamplesFromMAP(D0,D1,10000):')
	x = SamplesFromMAP(D0,D1,10000)
	print('mm = MarginalMomentsFromMAP(D0,D1,3):')
	mm = MarginalMomentsFromMAP(D0,D1,3)
	print('mm = ')
	print(mm)
	print("========================================")
	print('Testing BuTools function SamplesFromMMAP')
	print('Input:')
	print('------')
	D0 = ml.matrix([[-1.78, 0.29],[0.07, -0.92]])
	print('D0 = ')
	print(D0)
	D1 = ml.matrix([[0.15, 0.49],[0.23, 0.36]])
	print('D1 = ')
	print(D1)
	D2 = ml.matrix([[0.11, 0.2],[0.01, 0]])
	print('D2 = ')
	print(D2)
	D3 = ml.matrix([[0.14, 0.4],[0.11, 0.14]])
	print('D3 = ')
	print(D3)
	D = [D0,D1,D2,D3]
	print('Test:')
	print('-----')
	print('x = SamplesFromMMAP(D,10000):')
	x = SamplesFromMMAP(D,10000)
	print('mm = MarginalMomentsFromMMAP(D,3):')
	mm = MarginalMomentsFromMMAP(D,3)
	print('mm = ')
	print(mm)

if __name__ == "__main__":
	TestMAPPackage()

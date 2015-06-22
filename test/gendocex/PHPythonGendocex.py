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

def PHPythonGendocex():
	print('---BuTools: PH package test file---')
	print('Enable the verbose messages with the BuToolsVerbose flag')
	butools.verbose = True
	print('Enable input parameter checking with the BuToolsCheckInput flag')
	butools.checkInput = True
	np.set_printoptions(precision=5,linewidth=1024)
	outFile = open('PH_python.docex','w')
	with redirect_stdout(outFile):
		print('=== MomentsFromME ===')
		print('>>> a = ml.matrix([[0.2, 0.3, 0.5]])')
		a = ml.matrix([[0.2, 0.3, 0.5]])
		print('>>> A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])')
		A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])
		print('>>> moms = MomentsFromME(a,A)')
		moms = MomentsFromME(a,A)
		print('>>> print(moms)')
		print(moms)
		print('>>> moms = MomentsFromME(a,A,9)')
		moms = MomentsFromME(a,A,9)
		print('>>> print(moms)')
		print(moms)
		print('=== MomentsFromPH ===')
		print('>>> a = ml.matrix([[0.1, 0.9, 0]])')
		a = ml.matrix([[0.1, 0.9, 0]])
		print('>>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])')
		A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
		print('>>> moms = MomentsFromPH(a,A,5)')
		moms = MomentsFromPH(a,A,5)
		print('>>> print(moms)')
		print(moms)
		print('=== CdfFromME ===')
		print('>>> a = ml.matrix([[0.2, 0.3, 0.5]])')
		a = ml.matrix([[0.2, 0.3, 0.5]])
		print('>>> A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])')
		A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])
		print('>>> x = np.arange(0,5.01,0.01)')
		x = np.arange(0,5.01,0.01)
		print('>>> cdf = CdfFromME(a,A,x)')
		cdf = CdfFromME(a,A,x)
		print('>>> plt.plot(x,cdf)')
		print('=== CdfFromPH ===')
		print('>>> a = ml.matrix([[0.1, 0.9, 0]])')
		a = ml.matrix([[0.1, 0.9, 0]])
		print('>>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])')
		A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
		print('>>> x = np.arange(0,3.002,0.002)')
		x = np.arange(0,3.002,0.002)
		print('>>> cdf = CdfFromPH(a,A,x)')
		cdf = CdfFromPH(a,A,x)
		print('>>> plt.plot(x,cdf)')
		print('=== PdfFromME ===')
		print('>>> a = ml.matrix([[0.2, 0.3, 0.5]])')
		a = ml.matrix([[0.2, 0.3, 0.5]])
		print('>>> A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])')
		A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])
		print('>>> x = np.arange(0,5.01,0.01)')
		x = np.arange(0,5.01,0.01)
		print('>>> pdf = PdfFromME(a,A,x)')
		pdf = PdfFromME(a,A,x)
		print('>>> plt.plot(x,pdf)')
		print('=== PdfFromPH ===')
		print('>>> a = ml.matrix([[0.1, 0.9, 0]])')
		a = ml.matrix([[0.1, 0.9, 0]])
		print('>>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])')
		A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
		print('>>> x = np.arange(0,3.002,0.002)')
		x = np.arange(0,3.002,0.002)
		print('>>> pdf = PdfFromPH(a,A,x)')
		pdf = PdfFromPH(a,A,x)
		print('>>> plt.plot(x,pdf)')
		print('=== IntervalPdfFromPH ===')
		print('>>> a = ml.matrix([[0.1, 0.9, 0]])')
		a = ml.matrix([[0.1, 0.9, 0]])
		print('>>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])')
		A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
		print('>>> x = np.arange(0,3.002,0.002)')
		x = np.arange(0,3.002,0.002)
		print('>>> x,y = IntervalPdfFromPH(a,A,x)')
		x,y = IntervalPdfFromPH(a,A,x)
		print('>>> plt.plot(x,y)')
		print('=== RandomPH ===')
		print('>>> a,A = RandomPH(3,8,4)')
		a,A = RandomPH(3,8,4)
		print('>>> print(a)')
		print(a)
		print('>>> print(A)')
		print(A)
		print('=== CheckMERepresentation ===')
		print('>>> a = ml.matrix([[-0.2, 0.2]])')
		a = ml.matrix([[-0.2, 0.2]])
		print('>>> A = ml.matrix([[1, -1],[1, -2]])')
		A = ml.matrix([[1, -1],[1, -2]])
		print('>>> flag = CheckMERepresentation(a,A)')
		flag = CheckMERepresentation(a,A)
		print('>>> print(flag)')
		print(flag)
		print('>>> a = ml.matrix([[-0.2, 0.4, 0.8]])')
		a = ml.matrix([[-0.2, 0.4, 0.8]])
		print('>>> A = ml.matrix([[-2, 0, 3],[0, -1, 1],[0, -1, -1]])')
		A = ml.matrix([[-2, 0, 3],[0, -1, 1],[0, -1, -1]])
		print('>>> flag = CheckMERepresentation(a,A)')
		flag = CheckMERepresentation(a,A)
		print('>>> print(flag)')
		print(flag)
		print('>>> a = ml.matrix([[0.2, 0.3, 0.5]])')
		a = ml.matrix([[0.2, 0.3, 0.5]])
		print('>>> A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])')
		A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])
		print('>>> flag = CheckMERepresentation(a,A)')
		flag = CheckMERepresentation(a,A)
		print('>>> print(flag)')
		print(flag)
		print('=== CheckPHRepresentation ===')
		print('>>> a = ml.matrix([[0.2]])')
		a = ml.matrix([[0.2]])
		print('>>> A = ml.matrix([[-1, 1],[1, -2]])')
		A = ml.matrix([[-1, 1],[1, -2]])
		print('>>> flag = CheckPHRepresentation(a,A)')
		flag = CheckPHRepresentation(a,A)
		print('>>> print(flag)')
		print(flag)
		print('>>> a = ml.matrix([[0.2, 0.7]])')
		a = ml.matrix([[0.2, 0.7]])
		print('>>> A = ml.matrix([[-1, 1],[1, -2]])')
		A = ml.matrix([[-1, 1],[1, -2]])
		print('>>> flag = CheckPHRepresentation(a,A)')
		flag = CheckPHRepresentation(a,A)
		print('>>> print(flag)')
		print(flag)
		print('=== CheckMEPositiveDensity ===')
		print('>>> a = ml.matrix([[0.2, 0.3, 0.5]])')
		a = ml.matrix([[0.2, 0.3, 0.5]])
		print('>>> A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])')
		A = ml.matrix([[-1, 0, 0],[0, -3, 2],[0, -2, -3]])
		print('>>> flag = CheckMEPositiveDensity(a,A)')
		flag = CheckMEPositiveDensity(a,A)
		print('>>> print(flag)')
		print(flag)
		print('>>> a = ml.matrix([[0.2, 0.3, 0.5]])')
		a = ml.matrix([[0.2, 0.3, 0.5]])
		print('>>> A = ml.matrix([[-1, 0, 0],[0, -3, 2.9],[0, -2.9, -3]])')
		A = ml.matrix([[-1, 0, 0],[0, -3, 2.9],[0, -2.9, -3]])
		print('>>> flag = CheckMEPositiveDensity(a,A)')
		flag = CheckMEPositiveDensity(a,A)
		print('>>> print(flag)')
		print(flag)
		print('=== APHFrom3Moments ===')
		print('>>> moms = [10.0, 125.0, 8400.0]')
		moms = [10.0, 125.0, 8400.0]
		print('>>> a,A = APHFrom3Moments(moms)')
		a,A = APHFrom3Moments(moms)
		print('>>> print(a)')
		print(a)
		print('>>> print(A)')
		print(A)
		print('>>> phmoms = MomentsFromPH(a,A,3)')
		phmoms = MomentsFromPH(a,A,3)
		print('>>> print(phmoms)')
		print(phmoms)
		print('>>> moms = [10.0, 525.0, 31400.0]')
		moms = [10.0, 525.0, 31400.0]
		print('>>> a,A = APHFrom3Moments(moms)')
		a,A = APHFrom3Moments(moms)
		print('>>> print(a)')
		print(a)
		print('>>> print(A)')
		print(A)
		print('>>> phmoms = MomentsFromPH(a,A,3)')
		phmoms = MomentsFromPH(a,A,3)
		print('>>> print(phmoms)')
		print(phmoms)
		print('=== PH2From3Moments ===')
		print('>>> moms = [10.0, 160.0, 3500.0]')
		moms = [10.0, 160.0, 3500.0]
		print('>>> a,A = PH2From3Moments(moms)')
		a,A = PH2From3Moments(moms)
		print('>>> print(a)')
		print(a)
		print('>>> print(A)')
		print(A)
		print('>>> phmoms = MomentsFromPH(a,A,3)')
		phmoms = MomentsFromPH(a,A,3)
		print('>>> print(phmoms)')
		print(phmoms)
		print('>>> moms = [10.0, 260.0, 13500.0]')
		moms = [10.0, 260.0, 13500.0]
		print('>>> a,A = PH2From3Moments(moms)')
		a,A = PH2From3Moments(moms)
		print('>>> print(a)')
		print(a)
		print('>>> print(A)')
		print(A)
		print('>>> phmoms = MomentsFromPH(a,A,3)')
		phmoms = MomentsFromPH(a,A,3)
		print('>>> print(phmoms)')
		print(phmoms)
		print('=== PH3From5Moments ===')
		print('>>> a = ml.matrix([[0.1, 0.9, 0]])')
		a = ml.matrix([[0.1, 0.9, 0]])
		print('>>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])')
		A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
		print('>>> moms = MomentsFromPH(a,A)')
		moms = MomentsFromPH(a,A)
		print('>>> print(moms)')
		print(moms)
		print('>>> a,A = PH3From5Moments(moms)')
		a,A = PH3From5Moments(moms)
		print('>>> print(a)')
		print(a)
		print('>>> print(A)')
		print(A)
		print('>>> phmoms = MomentsFromME(a,A,5)')
		phmoms = MomentsFromME(a,A,5)
		print('>>> print(phmoms)')
		print(phmoms)
		print('>>> a = ml.matrix([[0.2, 0.3, 0.5]])')
		a = ml.matrix([[0.2, 0.3, 0.5]])
		print('>>> A = ml.matrix([[-1, 0, 0],[0, -3, 0.5],[0, -0.5, -3]])')
		A = ml.matrix([[-1, 0, 0],[0, -3, 0.5],[0, -0.5, -3]])
		print('>>> moms = MomentsFromME(a,A)')
		moms = MomentsFromME(a,A)
		print('>>> print(moms)')
		print(moms)
		print('>>> a,A = PH3From5Moments(moms)')
		a,A = PH3From5Moments(moms)
		print('>>> print(a)')
		print(a)
		print('>>> print(A)')
		print(A)
		print('>>> phmoms = MomentsFromME(a,A,5)')
		phmoms = MomentsFromME(a,A,5)
		print('>>> print(phmoms)')
		print(phmoms)
		print('=== MEFromMoments ===')
		print('>>> a = ml.matrix([[0.1, 0.9, 0]])')
		a = ml.matrix([[0.1, 0.9, 0]])
		print('>>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])')
		A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
		print('>>> moms = MomentsFromPH(a,A,5)')
		moms = MomentsFromPH(a,A,5)
		print('>>> print(moms)')
		print(moms)
		print('>>> a,A = MEFromMoments(moms)')
		a,A = MEFromMoments(moms)
		print('>>> print(a)')
		print(a)
		print('>>> print(A)')
		print(A)
		print('>>> memoms = MomentsFromME(a,A,5)')
		memoms = MomentsFromME(a,A,5)
		print('>>> print(memoms)')
		print(memoms)
		print('=== APH2ndMomentLowerBound ===')
		print('>>> mean = 1.9')
		mean = 1.9
		print('>>> n = 4')
		n = 4
		print('>>> mom2 = APH2ndMomentLowerBound(mean,n)')
		mom2 = APH2ndMomentLowerBound(mean,n)
		print('>>> print(mom2)')
		print(mom2)
		print('>>> cv2 = mom2/mean**2-1')
		cv2 = mom2/mean**2-1
		print('>>> print(1/cv2)')
		print(1/cv2)
		print('=== APH3rdMomentLowerBound ===')
		print('>>> mean = 1.9')
		mean = 1.9
		print('>>> mom2 = 5')
		mom2 = 5
		print('>>> n = 3')
		n = 3
		print('>>> mom3lower = APH3rdMomentLowerBound(mean,mom2,n)')
		mom3lower = APH3rdMomentLowerBound(mean,mom2,n)
		print('>>> print(mom3lower)')
		print(mom3lower)
		print('>>> mom3upper = APH3rdMomentUpperBound(mean,mom2,n)')
		mom3upper = APH3rdMomentUpperBound(mean,mom2,n)
		print('>>> print(mom3upper)')
		print(mom3upper)
		print('>>> mean = 1.9')
		mean = 1.9
		print('>>> mom2 = 5')
		mom2 = 5
		print('>>> n = 4')
		n = 4
		print('>>> mom3lower = APH3rdMomentLowerBound(mean,mom2,n)')
		mom3lower = APH3rdMomentLowerBound(mean,mom2,n)
		print('>>> print(mom3lower)')
		print(mom3lower)
		print('>>> mom3upper = APH3rdMomentUpperBound(mean,mom2,n)')
		mom3upper = APH3rdMomentUpperBound(mean,mom2,n)
		print('>>> print(mom3upper)')
		print(mom3upper)
		print('=== APH3rdMomentUpperBound ===')
		print('>>> mean = 1.9')
		mean = 1.9
		print('>>> mom2 = 5')
		mom2 = 5
		print('>>> n = 3')
		n = 3
		print('>>> mom3lower = APH3rdMomentLowerBound(mean,mom2,n)')
		mom3lower = APH3rdMomentLowerBound(mean,mom2,n)
		print('>>> print(mom3lower)')
		print(mom3lower)
		print('>>> mom3upper = APH3rdMomentUpperBound(mean,mom2,n)')
		mom3upper = APH3rdMomentUpperBound(mean,mom2,n)
		print('>>> print(mom3upper)')
		print(mom3upper)
		print('>>> mean = 1.9')
		mean = 1.9
		print('>>> mom2 = 5')
		mom2 = 5
		print('>>> n = 4')
		n = 4
		print('>>> mom3lower = APH3rdMomentLowerBound(mean,mom2,n)')
		mom3lower = APH3rdMomentLowerBound(mean,mom2,n)
		print('>>> print(mom3lower)')
		print(mom3lower)
		print('>>> mom3upper = APH3rdMomentUpperBound(mean,mom2,n)')
		mom3upper = APH3rdMomentUpperBound(mean,mom2,n)
		print('>>> print(mom3upper)')
		print(mom3upper)
		print('=== CanonicalFromPH2 ===')
		print('>>> a = ml.matrix([[0.12, 0.88]])')
		a = ml.matrix([[0.12, 0.88]])
		print('>>> A = ml.matrix([[-1.28, 0],[3.94, -3.94]])')
		A = ml.matrix([[-1.28, 0],[3.94, -3.94]])
		print('>>> b,B = CanonicalFromPH2(a,A)')
		b,B = CanonicalFromPH2(a,A)
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> Cm = SimilarityMatrix(A,B)')
		Cm = SimilarityMatrix(A,B)
		print('>>> err1 = la.norm(A*Cm-Cm*B)')
		err1 = la.norm(A*Cm-Cm*B)
		print('>>> err2 = la.norm(a*Cm-b)')
		err2 = la.norm(a*Cm-b)
		print('>>> print(np.max(err1,err2))')
		print(np.max(err1,err2))
		print('=== CanonicalFromPH3 ===')
		print('>>> a = ml.matrix([[0.1, 0.9, 0]])')
		a = ml.matrix([[0.1, 0.9, 0]])
		print('>>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])')
		A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
		print('>>> b,B = CanonicalFromPH3(a,A)')
		b,B = CanonicalFromPH3(a,A)
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> Cm = SimilarityMatrix(A,B)')
		Cm = SimilarityMatrix(A,B)
		print('>>> err1 = la.norm(A*Cm-Cm*B)')
		err1 = la.norm(A*Cm-Cm*B)
		print('>>> err2 = la.norm(a*Cm-b)')
		err2 = la.norm(a*Cm-b)
		print('>>> print(np.max(err1,err2))')
		print(np.max(err1,err2))
		print('=== AcyclicPHFromME ===')
		print('>>> a = ml.matrix([[-0.4, 1.4, 0.]])')
		a = ml.matrix([[-0.4, 1.4, 0.]])
		print('>>> A = ml.matrix([[-4., 1., 1.],[0., -2., 1.],[1., 0., -8.]])')
		A = ml.matrix([[-4., 1., 1.],[0., -2., 1.],[1., 0., -8.]])
		print('>>> b,B = AcyclicPHFromME(a,A)')
		b,B = AcyclicPHFromME(a,A)
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> ma = MomentsFromME(a,A,5)')
		ma = MomentsFromME(a,A,5)
		print('>>> print(ma)')
		print(ma)
		print('>>> mb = MomentsFromME(b,B,5)')
		mb = MomentsFromME(b,B,5)
		print('>>> print(mb)')
		print(mb)
		print('=== MonocyclicPHFromME ===')
		print('>>> a = ml.matrix([[0.2, 0.3, 0.5]])')
		a = ml.matrix([[0.2, 0.3, 0.5]])
		print('>>> A = ml.matrix([[-1., 0., 0.],[0., -3., 2.],[0., -2., -3.]])')
		A = ml.matrix([[-1., 0., 0.],[0., -3., 2.],[0., -2., -3.]])
		print('>>> b,B = MonocyclicPHFromME(a,A)')
		b,B = MonocyclicPHFromME(a,A)
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> ma = MomentsFromME(a,A,5)')
		ma = MomentsFromME(a,A,5)
		print('>>> print(ma)')
		print(ma)
		print('>>> mb = MomentsFromME(b,B,5)')
		mb = MomentsFromME(b,B,5)
		print('>>> print(mb)')
		print(mb)
		print('=== PHFromME ===')
		print('>>> a = ml.matrix([[-0.4, 1.4]])')
		a = ml.matrix([[-0.4, 1.4]])
		print('>>> A = ml.matrix([[-3.8, 2],[2, -9]])')
		A = ml.matrix([[-3.8, 2],[2, -9]])
		print('>>> flag = CheckMERepresentation(a,A)')
		flag = CheckMERepresentation(a,A)
		print('>>> print(flag)')
		print(flag)
		print('>>> flag = CheckPHRepresentation(a,A)')
		flag = CheckPHRepresentation(a,A)
		print('>>> print(flag)')
		print(flag)
		print('>>> b,B = PHFromME(a,A)')
		b,B = PHFromME(a,A)
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> flag = CheckPHRepresentation(b,B)')
		flag = CheckPHRepresentation(b,B)
		print('>>> print(flag)')
		print(flag)
		print('>>> Cm = SimilarityMatrix(A,B)')
		Cm = SimilarityMatrix(A,B)
		print('>>> err1 = la.norm(A*Cm-Cm*B)')
		err1 = la.norm(A*Cm-Cm*B)
		print('>>> err2 = la.norm(a*Cm-b)')
		err2 = la.norm(a*Cm-b)
		print('>>> print(np.max(err1,err2))')
		print(np.max(err1,err2))
		print('>>> a = ml.matrix([[-0.5, 1.5]])')
		a = ml.matrix([[-0.5, 1.5]])
		print('>>> A = ml.matrix([[-3.8, 2],[2, -9]])')
		A = ml.matrix([[-3.8, 2],[2, -9]])
		print('>>> flag = CheckMERepresentation(a,A)')
		flag = CheckMERepresentation(a,A)
		print('>>> print(flag)')
		print(flag)
		print('>>> flag = CheckPHRepresentation(a,A)')
		flag = CheckPHRepresentation(a,A)
		print('>>> print(flag)')
		print(flag)
		print('>>> b,B = PHFromME(a,A)')
		b,B = PHFromME(a,A)
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> flag = CheckPHRepresentation(b,B)')
		flag = CheckPHRepresentation(b,B)
		print('>>> print(flag)')
		print(flag)
		print('>>> Cm = SimilarityMatrix(A,B)')
		Cm = SimilarityMatrix(A,B)
		print('>>> err1 = la.norm(A*Cm-Cm*B)')
		err1 = la.norm(A*Cm-Cm*B)
		print('>>> err2 = la.norm(a*Cm-b)')
		err2 = la.norm(a*Cm-b)
		print('>>> print(np.max(err1,err2))')
		print(np.max(err1,err2))
		print('=== MEOrder ===')
		print('>>> a = ml.matrix([[1./6, 1./6, 1./6, 1./6, 1./6, 1./6]])')
		a = ml.matrix([[1./6, 1./6, 1./6, 1./6, 1./6, 1./6]])
		print('>>> A = ml.matrix([[-1., 0., 0., 0., 0., 0.],[0.5, -2., 1., 0., 0., 0.],[1., 0., -3., 1., 0., 0.],[1., 0., 1., -4., 1., 0.],[4., 0., 0., 0., -5., 0.],[5., 0., 0., 0., 0., -6.]])')
		A = ml.matrix([[-1., 0., 0., 0., 0., 0.],[0.5, -2., 1., 0., 0., 0.],[1., 0., -3., 1., 0., 0.],[1., 0., 1., -4., 1., 0.],[4., 0., 0., 0., -5., 0.],[5., 0., 0., 0., 0., -6.]])
		print('>>> co = MEOrder(a,A,"cont")')
		co = MEOrder(a,A,"cont")
		print('>>> print(co)')
		print(co)
		print('>>> oo = MEOrder(a,A,"obs")')
		oo = MEOrder(a,A,"obs")
		print('>>> print(oo)')
		print(oo)
		print('>>> coo = MEOrder(a,A,"obscont")')
		coo = MEOrder(a,A,"obscont")
		print('>>> print(coo)')
		print(coo)
		print('>>> mo = MEOrder(a,A,"moment")')
		mo = MEOrder(a,A,"moment")
		print('>>> print(mo)')
		print(mo)
		print('>>> a = ml.matrix([[2./3, 1./3]])')
		a = ml.matrix([[2./3, 1./3]])
		print('>>> A = ml.matrix([[-1., 1.],[0., -3.]])')
		A = ml.matrix([[-1., 1.],[0., -3.]])
		print('>>> co = MEOrder(a,A,"cont")')
		co = MEOrder(a,A,"cont")
		print('>>> print(co)')
		print(co)
		print('>>> oo = MEOrder(a,A,"obs")')
		oo = MEOrder(a,A,"obs")
		print('>>> print(oo)')
		print(oo)
		print('>>> coo = MEOrder(a,A,"obscont")')
		coo = MEOrder(a,A,"obscont")
		print('>>> print(coo)')
		print(coo)
		print('>>> mo = MEOrder(a,A,"moment")')
		mo = MEOrder(a,A,"moment")
		print('>>> print(mo)')
		print(mo)
		print('>>> b = ml.matrix([[0.2, 0.3, 0.5]])')
		b = ml.matrix([[0.2, 0.3, 0.5]])
		print('>>> B = ml.matrix([[-1., 0., 0.],[0., -3., 1.],[0., -1., -3.]])')
		B = ml.matrix([[-1., 0., 0.],[0., -3., 1.],[0., -1., -3.]])
		print('>>> a,A = MonocyclicPHFromME(b,B)')
		a,A = MonocyclicPHFromME(b,B)
		print('>>> print(a)')
		print(a)
		print('>>> print(A)')
		print(A)
		print('>>> co = MEOrder(a,A,"cont")')
		co = MEOrder(a,A,"cont")
		print('>>> print(co)')
		print(co)
		print('>>> oo = MEOrder(a,A,"obs")')
		oo = MEOrder(a,A,"obs")
		print('>>> print(oo)')
		print(oo)
		print('>>> coo = MEOrder(a,A,"obscont")')
		coo = MEOrder(a,A,"obscont")
		print('>>> print(coo)')
		print(coo)
		print('>>> mo = MEOrder(a,A,"moment")')
		mo = MEOrder(a,A,"moment")
		print('>>> print(mo)')
		print(mo)
		print('=== MEOrderFromMoments ===')
		print('>>> a = ml.matrix([[0.1, 0.9, 0]])')
		a = ml.matrix([[0.1, 0.9, 0]])
		print('>>> A = ml.matrix([[-6.2, 2., 0.],[2., -9., 1.],[1., 0., -3.]])')
		A = ml.matrix([[-6.2, 2., 0.],[2., -9., 1.],[1., 0., -3.]])
		print('>>> moms = MomentsFromME(a,A)')
		moms = MomentsFromME(a,A)
		print('>>> print(moms)')
		print(moms)
		print('>>> mo = MEOrderFromMoments(moms)')
		mo = MEOrderFromMoments(moms)
		print('>>> print(mo)')
		print(mo)
		print('>>> b = ml.matrix([[0.2, 0.3, 0.5]])')
		b = ml.matrix([[0.2, 0.3, 0.5]])
		print('>>> B = ml.matrix([[-1., 0., 0.],[0., -3., 2.],[0., -2., -3.]])')
		B = ml.matrix([[-1., 0., 0.],[0., -3., 2.],[0., -2., -3.]])
		print('>>> a,A = MonocyclicPHFromME(b,B)')
		a,A = MonocyclicPHFromME(b,B)
		print('>>> moms = MomentsFromME(a,A)')
		moms = MomentsFromME(a,A)
		print('>>> print(moms)')
		print(moms)
		print('>>> mo = MEOrderFromMoments(moms)')
		mo = MEOrderFromMoments(moms)
		print('>>> print(mo)')
		print(mo)
		print('=== MinimalRepFromME ===')
		print('>>> a = ml.matrix([[1./6, 1./6, 1./6, 1./6, 1./6, 1./6]])')
		a = ml.matrix([[1./6, 1./6, 1./6, 1./6, 1./6, 1./6]])
		print('>>> A = ml.matrix([[-1., 0., 0., 0., 0., 0.],[0.5, -2., 1., 0., 0., 0.],[1., 0., -3., 1., 0., 0.],[1., 0., 1., -4., 1., 0.],[4., 0., 0., 0., -5., 0.],[5., 0., 0., 0., 0., -6.]])')
		A = ml.matrix([[-1., 0., 0., 0., 0., 0.],[0.5, -2., 1., 0., 0., 0.],[1., 0., -3., 1., 0., 0.],[1., 0., 1., -4., 1., 0.],[4., 0., 0., 0., -5., 0.],[5., 0., 0., 0., 0., -6.]])
		print('>>> b,B = MinimalRepFromME(a,A,"cont")')
		b,B = MinimalRepFromME(a,A,"cont")
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> b,B = MinimalRepFromME(a,A,"obs")')
		b,B = MinimalRepFromME(a,A,"obs")
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> b,B = MinimalRepFromME(a,A,"obscont")')
		b,B = MinimalRepFromME(a,A,"obscont")
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> b,B = MinimalRepFromME(a,A,"moment")')
		b,B = MinimalRepFromME(a,A,"moment")
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> a = ml.matrix([[2./3, 1./3]])')
		a = ml.matrix([[2./3, 1./3]])
		print('>>> A = ml.matrix([[-1., 1.],[0., -3.]])')
		A = ml.matrix([[-1., 1.],[0., -3.]])
		print('>>> b,B = MinimalRepFromME(a,A,"cont")')
		b,B = MinimalRepFromME(a,A,"cont")
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> b,B = MinimalRepFromME(a,A,"obs")')
		b,B = MinimalRepFromME(a,A,"obs")
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> b,B = MinimalRepFromME(a,A,"obscont")')
		b,B = MinimalRepFromME(a,A,"obscont")
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> b,B = MinimalRepFromME(a,A,"moment")')
		b,B = MinimalRepFromME(a,A,"moment")
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> b = ml.matrix([[0.2, 0.3, 0.5]])')
		b = ml.matrix([[0.2, 0.3, 0.5]])
		print('>>> B = ml.matrix([[-1., 0., 0.],[0., -3., 1.],[0., -1., -3.]])')
		B = ml.matrix([[-1., 0., 0.],[0., -3., 1.],[0., -1., -3.]])
		print('>>> a,A = MonocyclicPHFromME(b,B)')
		a,A = MonocyclicPHFromME(b,B)
		print('>>> print(a)')
		print(a)
		print('>>> print(A)')
		print(A)
		print('>>> b,B = MinimalRepFromME(a,A,"cont")')
		b,B = MinimalRepFromME(a,A,"cont")
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> b,B = MinimalRepFromME(a,A,"obs")')
		b,B = MinimalRepFromME(a,A,"obs")
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> Cm = SimilarityMatrix(B,A)')
		Cm = SimilarityMatrix(B,A)
		print('>>> err1 = la.norm(B*Cm-Cm*A)')
		err1 = la.norm(B*Cm-Cm*A)
		print('>>> err2 = la.norm(b*Cm-a)')
		err2 = la.norm(b*Cm-a)
		print('>>> print(np.max(err1,err2))')
		print(np.max(err1,err2))
		print('>>> b,B = MinimalRepFromME(a,A,"obscont")')
		b,B = MinimalRepFromME(a,A,"obscont")
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> Cm = SimilarityMatrix(B,A)')
		Cm = SimilarityMatrix(B,A)
		print('>>> err1 = la.norm(B*Cm-Cm*A)')
		err1 = la.norm(B*Cm-Cm*A)
		print('>>> err2 = la.norm(b*Cm-a)')
		err2 = la.norm(b*Cm-a)
		print('>>> print(np.max(err1,err2))')
		print(np.max(err1,err2))
		print('>>> b,B = MinimalRepFromME(a,A,"moment")')
		b,B = MinimalRepFromME(a,A,"moment")
		print('>>> print(b)')
		print(b)
		print('>>> print(B)')
		print(B)
		print('>>> Cm = SimilarityMatrix(B,A)')
		Cm = SimilarityMatrix(B,A)
		print('>>> err1 = la.norm(B*Cm-Cm*A)')
		err1 = la.norm(B*Cm-Cm*A)
		print('>>> err2 = la.norm(b*Cm-a)')
		err2 = la.norm(b*Cm-a)
		print('>>> print(np.max(err1,err2))')
		print(np.max(err1,err2))
		print('=== SamplesFromPH ===')
		print('>>> a = ml.matrix([[0.1, 0.9, 0]])')
		a = ml.matrix([[0.1, 0.9, 0]])
		print('>>> A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])')
		A = ml.matrix([[-6.2, 2, 0],[2, -9, 1],[1, 0, -3]])
		print('>>> x = SamplesFromPH(a,A,1000)')
		x = SamplesFromPH(a,A,1000)
		print('>>> mp = MomentsFromPH(a,A,3)')
		mp = MomentsFromPH(a,A,3)
		print('>>> print(mp)')
		print(mp)
	outFile.close()

if __name__ == "__main__":
	PHPythonGendocex()

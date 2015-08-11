import sys
sys.path.append("/home/gabor/github/butools/Python")
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


print('---BuTools: Trace package test file---')
print('Enable the verbose messages with the BuToolsVerbose flag')
butools.verbose = True
print('Enable input parameter checking with the BuToolsCheckInput flag')
butools.checkInput = True
np.set_printoptions(precision=5,linewidth=1024)
print('========================================')
print('Testing BuTools function CdfFromTrace')
print('Input:')
print('------')
D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
print('D0 = ')
print(D0)
D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
print('D1 = ')
print(D1)
print('tr = SamplesFromMAP(D0, D1, 1000000):')
tr = SamplesFromMAP(D0, D1, 1000000)
print('Test:')
print('-----')
print('x, y = CdfFromTrace(tr):')
x, y = CdfFromTrace(tr)
plt.plot(x, y)
meandiff = np.abs(np.diff(x).dot(1.-y[0:-1])-np.mean(tr))/np.mean(tr)
assert np.all(np.diff(y)>=0) and y[-1]<=1 and y[0]>=0, "CdfFromTrace returned a wrong cdf!"
assert meandiff<10**-2, "The mean obtained from the cdf returned by CdfFromTrace does not match the trace mean!"
print('========================================')
print('Testing BuTools function PdfFromTrace')
print('Input:')
print('------')
D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
print('D0 = ')
print(D0)
D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
print('D1 = ')
print(D1)
x = np.arange(0.0,0.51,0.01)
print('tr = SamplesFromMAP(D0, D1, 1000000):')
tr = SamplesFromMAP(D0, D1, 1000000)
print('Test:')
print('-----')
print('x, y = PdfFromTrace(tr, x):')
x, y = PdfFromTrace(tr, x)
print('a, A = MarginalDistributionFromMAP(D0, D1):')
a, A = MarginalDistributionFromMAP(D0, D1)
print('xm, ym = IntervalPdfFromPH(a, A, x):')
xm, ym = IntervalPdfFromPH(a, A, x)
plt.plot(x, y, xm, ym)
meandiff = np.abs(x[0:-1].dot((np.diff(x)*y[0:-1]))-np.mean(tr))/np.mean(tr)
assert np.all(y>=0), "PdfFromTrace returned a wrong pdf!"
assert meandiff<10**-2, "The mean obtained from the pdf returned by PdfFromTrace does not match the trace mean!"
print('========================================')
print('Testing BuTools function MarginalMomentsFromTrace')
print('Input:')
print('------')
D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
print('D0 = ')
print(D0)
D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
print('D1 = ')
print(D1)
print('tr = SamplesFromMAP(D0, D1, 1000000):')
tr = SamplesFromMAP(D0, D1, 1000000)
print('Test:')
print('-----')
print('moms = MarginalMomentsFromTrace(tr, 3):')
moms = MarginalMomentsFromTrace(tr, 3)
print('moms = ')
print(moms)
print('mmoms = MarginalMomentsFromMAP(D0, D1, 3):')
mmoms = MarginalMomentsFromMAP(D0, D1, 3)
print('mmoms = ')
print(mmoms)
assert la.norm((np.array(moms)-np.array(mmoms))/np.array(mmoms))<10**-1, "Moments from MarginalMomentsFromTrace are far from the theoretical moments of the trace!"
print('========================================')
print('Testing BuTools function LagCorrelationsFromTrace')
print('Input:')
print('------')
D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
print('D0 = ')
print(D0)
D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
print('D1 = ')
print(D1)
print('tr = SamplesFromMAP(D0, D1, 1000000):')
tr = SamplesFromMAP(D0, D1, 1000000)
print('Test:')
print('-----')
print('acf = LagCorrelationsFromTrace(tr, 10):')
acf = LagCorrelationsFromTrace(tr, 10)
print('acf = ')
print(acf)
print('macf = LagCorrelationsFromMAP(D0, D1, 10):')
macf = LagCorrelationsFromMAP(D0, D1, 10)
print('macf = ')
print(macf)
assert la.norm(np.array(acf)-np.array(macf))<10**-1, "Autocorrelations from LagCorrelationsFromTrace are far from the theoretical autocorrelations of the trace!"
print('========================================')
print('Testing BuTools function LagkJointMomentsFromTrace')
print('Input:')
print('------')
D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
print('D0 = ')
print(D0)
D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
print('D1 = ')
print(D1)
print('tr = SamplesFromMAP(D0, D1, 1000000):')
tr = SamplesFromMAP(D0, D1, 1000000)
print('Test:')
print('-----')
print('Nm1 = LagkJointMomentsFromTrace(tr, 3, 1):')
Nm1 = LagkJointMomentsFromTrace(tr, 3, 1)
print('Nm1 = ')
print(Nm1)
print('mNm1 = LagkJointMomentsFromMAP(D0, D1, 3, 1):')
mNm1 = LagkJointMomentsFromMAP(D0, D1, 3, 1)
print('mNm1 = ')
print(mNm1)
assert la.norm((Nm1-mNm1)/mNm1)<10**-1, "Moments from LagkJointMomentsFromTrace are far from the theoretical moments of the trace!"
print('Test:')
print('-----')
print('Nm2 = LagkJointMomentsFromTrace(tr, 3, 2):')
Nm2 = LagkJointMomentsFromTrace(tr, 3, 2)
print('Nm2 = ')
print(Nm2)
print('mNm2 = LagkJointMomentsFromMAP(D0, D1, 3, 2):')
mNm2 = LagkJointMomentsFromMAP(D0, D1, 3, 2)
print('mNm2 = ')
print(mNm2)
assert la.norm((Nm2-mNm2)/mNm2)<10**-1, "Moments from LagkJointMomentsFromTrace are far from the theoretical moments of the trace!"
print('========================================')
print('Testing BuTools function CdfFromWeightedTrace')
print('Input:')
print('------')
wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34]
print('wtrace = ')
print(wtrace)
weights = [12., 1., 34., 23., 8., 2.]
print('weights = ')
print(weights)
print('Test:')
print('-----')
print('x, y = CdfFromWeightedTrace(wtrace, weights):')
x, y = CdfFromWeightedTrace(wtrace, weights)
plt.plot(x, y)
assert np.all(np.diff(y)>=0) and y[-1]<=1 and y[0]>=0, "CdfFromWeightedTrace returned a wrong cdf!"
print('========================================')
print('Testing BuTools function PdfFromWeightedTrace')
print('Input:')
print('------')
wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34]
print('wtrace = ')
print(wtrace)
weights = [12., 1., 34., 23., 8., 2.]
print('weights = ')
print(weights)
x = np.arange(0.0,3.1,0.1)
print('Test:')
print('-----')
print('x, y = PdfFromWeightedTrace(wtrace, weights, x):')
x, y = PdfFromWeightedTrace(wtrace, weights, x)
plt.plot(x, y)
assert np.all(y>=0), "PdfFromWeightedTrace returned a wrong pdf!"
print('========================================')
print('Testing BuTools function MarginalMomentsFromWeightedTrace')
print('Input:')
print('------')
wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34]
print('wtrace = ')
print(wtrace)
weights = [12., 1., 34., 23., 8., 2.]
print('weights = ')
print(weights)
print('Test:')
print('-----')
print('moms = MarginalMomentsFromWeightedTrace(wtrace, weights, 3):')
moms = MarginalMomentsFromWeightedTrace(wtrace, weights, 3)
print('moms = ')
print(moms)


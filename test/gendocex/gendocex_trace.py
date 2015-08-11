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
outFile = open('/home/gabor/github/butools/test/docex/Trace_python.docex','w')
with redirect_stdout(outFile):
    print('=== CdfFromTrace ===')
    print('>>> D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])')
    D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
    print('>>> D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])')
    D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
    print('>>> tr = SamplesFromMAP(D0, D1, 1000000)')
    tr = SamplesFromMAP(D0, D1, 1000000)
    print('>>> x, y = CdfFromTrace(tr)')
    x, y = CdfFromTrace(tr)
    print('>>> plt.plot(x, y)')
    meandiff = np.abs(np.diff(x).dot(1.-y[0:-1])-np.mean(tr))/np.mean(tr)
    print('=== PdfFromTrace ===')
    print('>>> D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])')
    D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
    print('>>> D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])')
    D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
    print('>>> x = np.arange(0.0,0.51,0.01)')
    x = np.arange(0.0,0.51,0.01)
    print('>>> tr = SamplesFromMAP(D0, D1, 1000000)')
    tr = SamplesFromMAP(D0, D1, 1000000)
    print('>>> x, y = PdfFromTrace(tr, x)')
    x, y = PdfFromTrace(tr, x)
    print('>>> a, A = MarginalDistributionFromMAP(D0, D1)')
    a, A = MarginalDistributionFromMAP(D0, D1)
    print('>>> xm, ym = IntervalPdfFromPH(a, A, x)')
    xm, ym = IntervalPdfFromPH(a, A, x)
    print('>>> plt.plot(x, y, xm, ym)')
    meandiff = np.abs(x[0:-1].dot((np.diff(x)*y[0:-1]))-np.mean(tr))/np.mean(tr)
    print('=== MarginalMomentsFromTrace ===')
    print('>>> D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])')
    D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
    print('>>> D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])')
    D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
    print('>>> tr = SamplesFromMAP(D0, D1, 1000000)')
    tr = SamplesFromMAP(D0, D1, 1000000)
    print('>>> moms = MarginalMomentsFromTrace(tr, 3)')
    moms = MarginalMomentsFromTrace(tr, 3)
    print('>>> print(moms)')
    print(moms)
    print('>>> mmoms = MarginalMomentsFromMAP(D0, D1, 3)')
    mmoms = MarginalMomentsFromMAP(D0, D1, 3)
    print('>>> print(mmoms)')
    print(mmoms)
    print('=== LagCorrelationsFromTrace ===')
    print('>>> D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])')
    D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
    print('>>> D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])')
    D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
    print('>>> tr = SamplesFromMAP(D0, D1, 1000000)')
    tr = SamplesFromMAP(D0, D1, 1000000)
    print('>>> acf = LagCorrelationsFromTrace(tr, 10)')
    acf = LagCorrelationsFromTrace(tr, 10)
    print('>>> print(acf)')
    print(acf)
    print('>>> macf = LagCorrelationsFromMAP(D0, D1, 10)')
    macf = LagCorrelationsFromMAP(D0, D1, 10)
    print('>>> print(macf)')
    print(macf)
    print('=== LagkJointMomentsFromTrace ===')
    print('>>> D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])')
    D0 = ml.matrix([[-18., 1., 4.],[2., -18., 7.],[1., 3., -32.]])
    print('>>> D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])')
    D1 = ml.matrix([[12., 1., 0.],[1., 8., 0.],[2., 1., 25.]])
    print('>>> tr = SamplesFromMAP(D0, D1, 1000000)')
    tr = SamplesFromMAP(D0, D1, 1000000)
    print('>>> Nm1 = LagkJointMomentsFromTrace(tr, 3, 1)')
    Nm1 = LagkJointMomentsFromTrace(tr, 3, 1)
    print('>>> print(Nm1)')
    print(Nm1)
    print('>>> mNm1 = LagkJointMomentsFromMAP(D0, D1, 3, 1)')
    mNm1 = LagkJointMomentsFromMAP(D0, D1, 3, 1)
    print('>>> print(mNm1)')
    print(mNm1)
    print('>>> Nm2 = LagkJointMomentsFromTrace(tr, 3, 2)')
    Nm2 = LagkJointMomentsFromTrace(tr, 3, 2)
    print('>>> print(Nm2)')
    print(Nm2)
    print('>>> mNm2 = LagkJointMomentsFromMAP(D0, D1, 3, 2)')
    mNm2 = LagkJointMomentsFromMAP(D0, D1, 3, 2)
    print('>>> print(mNm2)')
    print(mNm2)
    print('=== CdfFromWeightedTrace ===')
    print('>>> wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34]')
    wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34]
    print('>>> weights = [12., 1., 34., 23., 8., 2.]')
    weights = [12., 1., 34., 23., 8., 2.]
    print('>>> x, y = CdfFromWeightedTrace(wtrace, weights)')
    x, y = CdfFromWeightedTrace(wtrace, weights)
    print('>>> plt.plot(x, y)')
    print('=== PdfFromWeightedTrace ===')
    print('>>> wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34]')
    wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34]
    print('>>> weights = [12., 1., 34., 23., 8., 2.]')
    weights = [12., 1., 34., 23., 8., 2.]
    print('>>> x = np.arange(0.0,3.1,0.1)')
    x = np.arange(0.0,3.1,0.1)
    print('>>> x, y = PdfFromWeightedTrace(wtrace, weights, x)')
    x, y = PdfFromWeightedTrace(wtrace, weights, x)
    print('>>> plt.plot(x, y)')
    print('=== MarginalMomentsFromWeightedTrace ===')
    print('>>> wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34]')
    wtrace = [0.12, 1.23, 0.546, 0.6765, 1.34, 2.34]
    print('>>> weights = [12., 1., 34., 23., 8., 2.]')
    weights = [12., 1., 34., 23., 8., 2.]
    print('>>> moms = MarginalMomentsFromWeightedTrace(wtrace, weights, 3)')
    moms = MarginalMomentsFromWeightedTrace(wtrace, weights, 3)
    print('>>> print(moms)')
    print(moms)


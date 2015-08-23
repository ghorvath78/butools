butools.fitting.EmpiricalRelativeEntropy
========================================

.. currentmodule:: butools.fitting

.. np:function:: EmpiricalRelativeEntropy

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`re = EmpiricalRelativeEntropy(f1, f2, intBounds)`
        * - Mathematica:
          - :code:`re = EmpiricalRelativeEntropy[f1, f2, intBounds]`
        * - Python/Numpy:
          - :code:`re = EmpiricalRelativeEntropy(f1, f2, intBounds)`


    Returns the relative entropy (aka Kullback–Leibler
    divergence) of two continuous functions given by samples
    and the bounds of the corresponding intervalls.
    
    This function can be used to characterize the distance
    between two density functions, distribution functions, 
    etc.
    
    Parameters
    ----------
    f1 : vector, length M
        Samples of the first continuous function
    f2 : vector, length M
        Samples of the second continuous function
    intBounds : vector, length M+1
        The bounds of the intervals. The ith sample
        corresponds to the interval 
        (intbounds(i),intbounds(i+1))

    Returns
    -------
    re : double
        The relative entropy

    Examples
    ========
    For Matlab:

    >>> tr = dlmread('/home/gabor/github/butools/test/data/bctrace.iat');
    >>> intBounds = linspace(0, MarginalMomentsFromTrace(tr, 1)*4, 50);
    >>> [pdfTrX, pdfTrY] = PdfFromTrace(tr, intBounds);
    >>> [cdfTrX, cdfTrY] = CdfFromTrace(tr);
    >>> step = ceil(length(tr)/2000);
    >>> cdfTrX = cdfTrX(1:step:length(tr));
    >>> cdfTrY = cdfTrY(1:step:length(tr));
    >>> [alpha, A] = APHFrom3Moments(MarginalMomentsFromTrace(tr, 3));
    >>> [pdfPHX, pdfPHY] = IntervalPdfFromPH(alpha, A, intBounds);
    >>> cdfPHY = CdfFromPH(alpha, A, cdfTrX);
    >>> rePdf = EmpiricalRelativeEntropy(pdfTrY, pdfPHY, intBounds);
    >>> disp(rePdf);
          0.43241
    >>> reCdf = EmpiricalRelativeEntropy(cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX);
    >>> disp(reCdf);
       0.00040609

    For Mathematica:

    
    For Python/Numpy:

    >>> tr = np.loadtxt("/home/gabor/github/butools/test/data/bctrace.iat")
    >>> intBounds = np.linspace(0, MarginalMomentsFromTrace(tr, 1)[0]*4, 50)
    >>> pdfTrX, pdfTrY = PdfFromTrace(tr, intBounds)
    >>> cdfTrX, cdfTrY = CdfFromTrace(tr)
    >>> step = math.ceil(Length(tr)/2000)
    >>> cdfTrX = cdfTrX[0:Length(tr):step]
    >>> cdfTrY = cdfTrY[0:Length(tr):step]
    >>> alpha, A = APHFrom3Moments(MarginalMomentsFromTrace(tr, 3))
    >>> pdfPHX, pdfPHY = IntervalPdfFromPH(alpha, A, intBounds)
    >>> cdfPHY = CdfFromPH(alpha, A, cdfTrX)
    >>> rePdf = EmpiricalRelativeEntropy(pdfTrY, pdfPHY, intBounds)
    >>> print(rePdf)
    0.432414379777
    >>> reCdf = EmpiricalRelativeEntropy(cdfTrY[0:-1], cdfPHY[0:-1], cdfTrX)
    >>> print(reCdf)
    0.000406094874315


butools.fitting.EmpiricalSquaredDifference
==========================================

.. currentmodule:: butools.fitting

.. np:function:: EmpiricalSquaredDifference

    .. list-table::
        :widths: 25 150

        * - Matlab:
          - :code:`sd = EmpiricalSquaredDifference(f1, f2, intBounds)`
        * - Mathematica:
          - :code:`sd = EmpiricalSquaredDifference[f1, f2, intBounds]`
        * - Python/Numpy:
          - :code:`sd = EmpiricalSquaredDifference(f1, f2, intBounds)`


    Returns the squared difference of two continuous 
    functions given by samples and the bounds of the 
    corresponding intervalls.
    
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
    sd : double
        The squared difference

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
    >>> sqPdf = EmpiricalSquaredDifference(pdfTrY, pdfPHY, intBounds);
    >>> disp(sqPdf);
         0.011854
    >>> sqCdf = EmpiricalSquaredDifference(cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX);
    >>> disp(sqCdf);
       3.8247e-10

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
    >>> sqPdf = EmpiricalSquaredDifference(pdfTrY, pdfPHY, intBounds)
    >>> print(sqPdf)
    0.0118541986064
    >>> sqCdf = EmpiricalSquaredDifference(cdfTrY[0:-1], cdfPHY[0:-1], cdfTrX)
    >>> print(sqCdf)
    3.8246917213e-10


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
    --------    
    For Matlab:
    
    >>> tr = dlmread('trace.txt');
    >>> intBounds = linspace(0, MarginalMomentsFromTrace(tr,1)*4, 50);
    >>> [pdfTrX, pdfTrY] = PdfFromTrace(tr,intBounds);
    >>> [alpha,A]=PHFromTrace(tr,5)   
    >>> [pdfPHX, pdfPHY] = IntervalPdfFromPH(alpha, A, intBounds);
    >>> EmpiricalRelativeEntropy (pdfTrY, pdfPHY, intBounds)
          0.37827
    

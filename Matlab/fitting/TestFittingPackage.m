function TestFittingPackage

    format short g;

    disp('---BuTools: Continous PH function test file---');

    disp('Enable the verbose messages with the BuToolsVerbose flag')
    global BuToolsVerbose;
    BuToolsVerbose = true

    disp('Enable input parameter checking with the BuToolsCheckInput flag')
    global BuToolsCheckInput;
    BuToolsCheckInput = true

    unzip('bctrace.zip');
    tr = dlmread('bctrace.iat');
    
    fprintf('Length of the trace: %d samples\n', length(tr));

disp('----------------------------------------------------------------------------');
    
    disp('Generating input for the test');
    disp('=============================');
    
    intBounds = linspace(0, MarginalMomentsFromTrace(tr,1)*4, 50);
    disp('Obtaining pdf of the trace: [pdfTrX, pdfTrY] = PdfFromTrace(tr,intBounds)');
    [pdfTrX, pdfTrY] = PdfFromTrace(tr,intBounds);
    disp('Obtaining cdf of the trace: [cdfTrX, cdfTrY] = CdfFromTrace(tr)');
    [cdfTrX, cdfTrY] = CdfFromTrace(tr);
    
    maxCdfPoints = 2000;
    if length(tr)>maxCdfPoints
        step = ceil(length(tr) / maxCdfPoints);
        ix = 1:step:length(tr);
        cdfTrX = cdfTrX(ix);
        cdfTrY = cdfTrY(ix);
    end
    
    disp('Fitting an APH based on 3 moments to the trace...');
    moms = MarginalMomentsFromTrace(tr,3);
    [alpha,A] = APHFrom3Moments(moms)
    disp('Obtaining pdf of the APH(3): [pdfPHX, pdfPHY] = IntervalPdfFromPH(alpha, A, intBounds)');
    [pdfPHX, pdfPHY] = IntervalPdfFromPH(alpha, A, intBounds);
    disp('Obtaining cdf of the APH(3): cdfPHY = CdfFromPH(alpha, A, cdfTrX)');
    cdfPHY = CdfFromPH(alpha, A, cdfTrX);
    
disp('----------------------------------------------------------------------------');

    disp('Testing distance functions');
    disp('==========================');
    
    disp('Calculating pdf squared difference: sqPdf = EmpiricalSquaredDifference (pdfTrY, pdfPHY, intBounds)');
    sqPdf = EmpiricalSquaredDifference (pdfTrY, pdfPHY, intBounds)
    
    disp('Calculating pdf relative entropy: rePdf = EmpiricalRelativeEntropy (pdfTrY, pdfPHY, intBounds)');
    rePdf = EmpiricalRelativeEntropy (pdfTrY, pdfPHY, intBounds)
    
    disp('Calculating cdf squared difference: sqCdf = EmpiricalSquaredDifference (cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX)');
    sqCdf = EmpiricalSquaredDifference (cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX)
    
    disp('Calculating cdf relative entropy: reCdf = EmpiricalRelativeEntropy (cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX)');
    reCdf = EmpiricalRelativeEntropy (cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX)
    
    disp('Calculating likelihood: logli = LikelihoodFromTrace(tr,alpha,A)');   
    logli = LikelihoodFromTrace(tr,alpha,A)
    
disp('----------------------------------------------------------------------------');

    disp('Fitting a PH(5) by using G-FIT');
    disp('==============================');
    
    disp('Perform fitting: [alpha,A]=PHFromTrace(tr,5)');   
    [alpha,A]=PHFromTrace(tr,5)
    
    [pdfPHX, pdfPHY] = IntervalPdfFromPH(alpha, A, intBounds, 1e-12);
    cdfPHY = CdfFromPH(alpha, A, cdfTrX, 1e-12);

    disp('Calculating pdf squared difference: sqPdf = EmpiricalSquaredDifference (pdfTrY, pdfPHY, intBounds)');
    sqPdf = EmpiricalSquaredDifference (pdfTrY, pdfPHY, intBounds)
    
    disp('Calculating pdf relative entropy: rePdf = EmpiricalRelativeEntropy (pdfTrY, pdfPHY, intBounds)');
    rePdf = EmpiricalRelativeEntropy (pdfTrY, pdfPHY, intBounds)
    
    disp('Calculating cdf squared difference: sqCdf = EmpiricalSquaredDifference (cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX)');
    sqCdf = EmpiricalSquaredDifference (cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX)
    
    disp('Calculating cdf relative entropy: reCdf = EmpiricalRelativeEntropy (cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX)');
    reCdf = EmpiricalRelativeEntropy (cdfTrY(1:end-1), cdfPHY(1:end-1), cdfTrX)
    
    disp('Calculating likelihood: logli = LikelihoodFromTrace(tr,alpha,A)');   
    logli = LikelihoodFromTrace(tr,alpha,A)
    
disp('----------------------------------------------------------------------------');

    disp('Testing further distance functions');
    disp('==================================');

    maxMAPTraceLen = 10000;
    if length(tr)>maxMAPTraceLen
        tr = tr(1:maxMAPTraceLen);
    end

    disp('Fitting a MAP based on 3 moments and 1 lag correlation...');
    corr1 = LagCorrelationsFromTrace(tr,1);
    [D0,D1]=MAPFromFewMomentsAndCorrelations(moms,corr1)

    disp('Calculating likelihood: logli = LikelihoodFromTrace(tr,D0,D1)');   
    logli = LikelihoodFromTrace(tr,D0,D1)
    
    disp('Obtaining acf of the trace: trAcf = LagCorrelationsFromTrace(tr, 10)');
    trAcf = LagCorrelationsFromTrace(tr, 10)
    
    disp('Obtaining acf of the MAP: mAcf = LagCorrelationsFromMAP(D0, D1, 10)');
    mAcf = LagCorrelationsFromMAP(D0, D1, 10, 1e-13)

    disp('Calculating acf squared difference: sqAcf = SquaredDifference (mAcf, trAcf)');
    sqAcf = SquaredDifference (mAcf, trAcf)
    
    disp('Calculating acf relative entropy: reAcf = RelativeEntropy (mAcf, trAcf)');
    reAcf = RelativeEntropy (mAcf, trAcf)

    
    disp('Fitting a MAP(5) by using SPEM-FIT');
    disp('==================================');
    
    disp('Perform fitting: [D0,D1]=MAPFromTrace(tr,5)');
    [D0,D1]=MAPFromTrace(tr,5)
    
    disp('Calculating likelihood: logli = LikelihoodFromTrace(tr,D0,D1)');   
    logli = LikelihoodFromTrace(tr,D0,D1)

    disp('Obtaining acf of the MAP: mAcf = LagCorrelationsFromMAP(D0, D1, 10)');
    mAcf = LagCorrelationsFromMAP(D0, D1, 10, 1e-13)

    disp('Calculating acf squared difference: sqAcf = SquaredDifference (mAcf, trAcf)');
    sqAcf = SquaredDifference (mAcf, trAcf)
    
    disp('Calculating acf relative entropy: reAcf = RelativeEntropy (mAcf, trAcf)');
    reAcf = RelativeEntropy (mAcf, trAcf)
end


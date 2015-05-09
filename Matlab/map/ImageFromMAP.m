%  ImageFromMAP(D0, D1, outFileName, prec)
%  
%  Depicts the given Markovian arrival process, and either
%  displays it or saves it to file.
%  
%  Parameters
%  ----------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the Markovian arrival process
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the Markovian arrival process
%  outFileName : string, optional
%      If it is not provided, or equals to 'display', the
%      image is displayed on the screen, otherwise it is 
%      written to the file. The file format is deduced 
%      from the file name.
%  prec : double, optional
%      Transition rates less then prec are considered to
%      be zero and are left out from the image. The 
%      default value is 1e-13.
%  
%  Notes
%  -----
%  The 'graphviz' software must be installed and available
%  in the path to use this feature.

function ImageFromMAP(D0,D1,outFileName,prec)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMAPRepresentation(D0,D1)
        error('ImageFromMAP: Input isn''t a valid MAP representation!');
    end
    
    if ~exist('prec','var')
        prec = 1e-13;
    end

    if ~exist('outFileName','var')
        outFileName = 'display';
    end
    
    ImageFromMMAP({D0,D1},outFileName,prec);
end


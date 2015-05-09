%  D = MMAPFromMRAP(H, precision)
%  
%  Obtains a Markovian representation of a rational
%  arrival process of the same size, if possible, using the
%  procedure published in [1]_.
%  
%  Parameters
%  ----------
%  H : list/cell of matrices of shape(M,M), length(N)
%      The H0...HN matrices of the MRAP to transform
%  precision : double, optional
%      A representation is considered to be a Markovian one
%      if it is closer to it than this precision
%  
%  Returns
%  -------
%  D : list/cell of matrices of shape(M,M), length(N)
%      The D0...DN matrices of the MMAP (if found)
%  
%  References
%  ----------
%  .. [1] András Horváth, Gábor Horváth, Miklós Telek, "A 
%         traffic based decomposition of two-class queueing 
%         networks with priority service". COMPUTER NETWORKS 
%         53:(8) pp. 1235-1248. (2009)

function D = MMAPFromMRAP (H,prec)

    if ~exist('prec','var')
        prec = 1e-14;
    end

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMRAPRepresentation(H)
        error('MMAPFromMRAP: Input isn''t a valid MRAP representation');
    end

    function nH = transfun (oH, B)
        nH = cell(1,length(oH));
        for i=1:length(oH)
            nH{i} = inv(B)*oH{i}*B;
        end
    end
        
    function dist = evalfun (oH, k)
        if nargin<2
            k = 0;
        end
        oH0 = oH{1} - diag(diag(oH{1}));
        if rem(k,2) == 0
            dist = min(min(oH0));
            for k=2:length(oH)
                dist = min(dist, min(min(oH{k})));
            end
        else
            dist = sum(oH0(oH0<0));
            for k=2:length(oH)
                oHk = oH{k};
                dist = dist + sum(oHk(oHk<0));
            end
        end
        dist = -dist;
    end

    D = FindMarkovianRepresentation (H, @transfun, @evalfun, prec);
end

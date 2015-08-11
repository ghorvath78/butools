%  r = CheckMoments(m, prec)
%  
%  Checks if the given moment sequence is valid in the sense
%  that it belongs to a distribution with support (0,inf).
%  
%  This procedure checks the determinant of `\Delta_n`
%  and `\Delta_n^{(1)}` according to [1]_.
%  
%  Parameters
%  ----------
%  m : list of doubles, length 2N+1
%      The (raw) moments to check 
%      (starts with the first moment).
%      Its length must be odd.
%  prec : double, optional
%      Entries with absolute value less than prec are 
%      considered to be zeros. The default value is 1e-14.
%      
%  Returns
%  -------
%  r : bool
%      The result of the check.
%  
%  References
%  ----------
%  .. [1] http://en.wikipedia.org/wiki/Stieltjes_moment_problem

function r = CheckMoments (m, prec)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end
    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-14;
    end
    
    if ~exist('prec','var')
        prec=BuToolsCheckPrecision;
    end

    if BuToolsCheckInput && mod(length(m),2)==0
        error ('CheckMoments: the number of moments must be odd!');
    end
    
    m = reshape(m, 1, length(m));
    m = [1.0 m];
    N = floor(length(m)/2)-1;
    
    for n=0:N
        H = hankel(m(1:n+1), m(n+1:2*n+1));
        H0 = hankel(m(2:n+2), m(n+2:2*n+2));
        if det(H)<-prec || det(H0)<-prec
            r = false;
            return;
        end
    end
    r = true;
end


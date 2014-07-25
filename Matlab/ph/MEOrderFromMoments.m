%  order = MEOrderFromMoments(moms, prec)
%  
%  Returns the order of ME distribution that can realize
%  the given moments.
%  
%  Parameters
%  ----------
%  moms : list of doubles
%      The list of moments
%  prec : double, optional
%      Precision used to detect if the determinant of the
%      Hankel matrix is zero. The default value is 1e-12.
%  
%  Returns
%  -------
%  order : int
%      The order of ME distribution that can realize the 
%      given moments
%  
%  References
%  ----------
%  .. [1]  L. Bodrog, A. Horvath, M. Telek, "Moment 
%          characterization of matrix exponential and Markovian
%          arrival processes," Annals of Operations Research, 
%          vol. 160, pp. 51-68, 2008.

function order = MEOrderFromMoments (moms, prec)

    if ~exist('prec','var')
        prec = 1e-12;
    end
    
    sizem = floor((length(moms)+1)/2);
    rmoms = [1 ReducedMomsFromMoms(moms)];
    
    for k=1:sizem
        hankel = zeros(k);
        for i=1:k
           for j=1:k
               hankel(i,j) = rmoms(i+j-1);
           end
        end
        if abs(det(hankel)) < prec
           order = k-1;
           return;
        end
    end
    order = sizem;
end

%  hm = HankelMomsFromMoms(m)
%  
%  Returns the Hankel moments given the raw moments.
%  
%  The raw moments are: `m_i=E(\mathcal{X}^i)`
%  
%  The ith Hankel moment is the determinant of matrix 
%  `\Delta_{i/2}`, if i is even, 
%  and it is the determinant of `\Delta^{(1)}_{(i+1)/2}`, 
%  if i is odd. For the definition of matrices `\Delta`
%  and `\Delta^{(1)}` see [1]_.
%     
%  Parameters
%  ----------
%  m : vector of doubles
%      The list of raw moments (starting with the first
%      moment)
%      
%  Returns
%  -------
%  hm : vector of doubles
%      The list of Hankel moments
%  
%  References
%  ----------
%  .. [1] http://en.wikipedia.org/wiki/Stieltjes_moment_problem

function hm=HankelMomsFromMoms(m)

    hm = [];
    for i=0:length(m)-1
        if rem(i,2) == 0
            N = i/2 + 1;
            H = hankel(m(1:N), m(N:2*N-1));
        else
            N = (i+1) / 2 + 1;
            H = hankel([1,m(1:N-1)], m(N-1:2*N-2));
        end
        hm = [hm, det(H)];
    end
end
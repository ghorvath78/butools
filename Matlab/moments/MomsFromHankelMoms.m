%  m = MomsFromHankelMoms(hm)
%  
%  Returns the raw moments given the Hankel moments.
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
%  hm : vector of doubles
%      The list of Hankel moments (starting with the first
%      moment)
%      
%  Returns
%  -------
%  m : vector of doubles
%      The list of raw moments
%  
%  References
%  ----------
%  .. [1] http://en.wikipedia.org/wiki/Stieltjes_moment_problem

function m=MomsFromHankelMoms(hm)

    m = hm(1);
    for i=1:length(hm)-1
        if rem(i,2) == 0
            N = i/2 + 1;
            H = hankel(m(1:N), [m(N:2*N-2),0]);
        else
            N = (i+1) / 2 + 1;
            H = hankel([1,m(1:N-1)], [m(N-1:2*N-3),0]);
        end
        h = hm(i+1);
        rH = H(1:N-1,:);
        for j=0:N-1
            rHd = rH;
            rHd(:,j+1) = [];
            cofactor = (-1)^(N+j-1) * det(rHd);
            if j<N-1
                h = h - cofactor * H(N,j+1);
            else
                m = [m h/cofactor];
            end
        end
    end
end


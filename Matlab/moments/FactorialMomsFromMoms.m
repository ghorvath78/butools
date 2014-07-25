%  fm = FactorialMomsFromMoms(m)
%  
%  Returns the factorial moments given the raw moments.
%  
%  The raw moments are: `m_i=E(\mathcal{X}^i)`
%  
%  The factorial moments are: `f_i=E(\mathcal{X}(\mathcal{X}-1)\cdots(\mathcal{X}-i+1))`
%     
%  Parameters
%  ----------
%  m : vector of doubles
%      The list of raw moments (starting with the first
%      moment)
%      
%  Returns
%  -------
%  fm : vector of doubles
%      The list of factorial moments
%  
%  References
%  ----------
%  http://en.wikipedia.org/wiki/Factorial_moment    

function fm=FactorialMomsFromMoms(m)

    n=length(m);
    fm=zeros(size(m));
    m = reshape(m, 1, n);
    for i=1:n
        eh=poly(0:i-1);
        eh=eh(end-1:-1:1);
        fm(i)=eh*(m(1:i))';    
    end
end

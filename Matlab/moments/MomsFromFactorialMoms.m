%  m = MomsFromFactorialMoms(fm)
%  
%  Returns the raw moments given the factorial moments.
%  
%  The raw moments are: `m_i=E(\mathcal{X}^i)`
%  
%  The factorial moments are: `f_i=E(\mathcal{X}(\mathcal{X}-1)\cdots(\mathcal{X}-i+1))`
%     
%  Parameters
%  ----------
%  fm : vector of doubles
%      The list of factorial moments (starting with the first
%      moment)
%      
%  Returns
%  -------
%  m : vector of doubles
%      The list of raw moments
%  
%  References
%  ----------
%  http://en.wikipedia.org/wiki/Factorial_moment    

function m=MomsFromFactorialMoms(fm)

    n=length(fm);
    m=zeros(1,n);
    m(1)=fm(1);

    for i=2:n
        eh=-poly(0:i-1);
        eh=eh(end-1:-1:2);
        m(i)=fm(i)+eh*(m(1:i-1))';
    end
    m = reshape(m, size(fm));
end

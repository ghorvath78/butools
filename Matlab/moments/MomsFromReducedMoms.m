%  m = MomsFromReducedMoms(rm)
%  
%  Returns the raw moments given the reduced moments.
%  
%  The raw moments are: `m_i=E(\mathcal{X}^i)`
%  
%  The reduced moments are: `\displaystyle r_i=\frac{m_i}{i!}`
%     
%  Parameters
%  ----------
%  rm : vector of doubles
%      The list of reduced moments (starting with the first
%      moment)
%      
%  Returns
%  -------
%  m : vector of doubles
%      The list of raw moments

function m=MomsFromReducedMoms(rm)

    m = rm;
    f = 1.0;
    for i=1:length(m)
        f = f * i;
        m(i) = m(i) * f;
    end
end
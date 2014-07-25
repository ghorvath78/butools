%  rm = ReducedMomsFromMoms(m)
%  
%  Returns the reduced moments given the raw moments.
%  
%  The raw moments are: `m_i=E(\mathcal{X}^i)`
%  
%  The reduced moments are: `\displaystyle r_i=\frac{m_i}{i!}`
%     
%  Parameters
%  ----------
%  m : vector of doubles
%      The list of raw moments (starting with the first
%      moment)
%      
%  Returns
%  -------
%  rm : vector of doubles
%      The list of reduced moments

function rm=ReducedMomsFromMoms(m)

    rm = m;
    f = 1.0;
    for i=1:length(m)
        f = f / i;
        rm(i) = rm(i) * f;
    end
end

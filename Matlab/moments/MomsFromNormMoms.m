%  m = MomsFromNormMoms(nm)
%  
%  Returns the raw moments given the normalized moments.
%  
%  The raw moments are: `m_i=E(\mathcal{X}^i)`
%  
%  The normalized moments are: `\displaystyle n_i=\frac{m_i}{m_{i-1} m_1}`
%     
%  Parameters
%  ----------
%  nm : vector of doubles
%      The list of normalized moments (starting with the first
%      moment)
%      
%  Returns
%  -------
%  m : vector of doubles
%      The list of raw moments

function m=MomsFromNormMoms(nm)

    m = nm;
    for i=2:length(nm)
        m(i) = m(1) * m(i) * m(i-1);
    end
end

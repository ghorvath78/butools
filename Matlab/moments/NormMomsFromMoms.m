%  nm = NormMomsFromMoms(m)
%  
%  Returns the normalized moments given the raw moments.
%  
%  The raw moments are: `m_i=E(\mathcal{X}^i)`
%  
%  The normalized moments are: `\displaystyle n_i=\frac{m_i}{m_{i-1} m_1}`
%     
%  Parameters
%  ----------
%  m : vector of doubles
%      The list of raw moments (starting with the first
%      moment)
%      
%  Returns
%  -------
%  nm : vector of doubles
%      The list of normalized moments

function nm=NormMomsFromMoms(m)

    nm = zeros(size(m));
    nm(1) = m(1);
    for i=2:length(m)
        nm(i) = m(i) / m(i-1) / m(1);
    end
end
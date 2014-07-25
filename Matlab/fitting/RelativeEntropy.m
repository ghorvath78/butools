%  re = RelativeEntropy(p1, p2)
%  
%  Returns the relative entropy (aka Kullback–Leibler
%  divergence) of two vectors.
%  
%  Parameters
%  ----------
%  p1 : vector, length M
%      The first vector
%  p2 : vector, length M
%      The second vector
%  
%  Returns
%  -------
%  re : double
%      The relative entropy calculated as
%      `re=\sum_{i=1}^M p1_i |\log(p1_i/p2_i)|`

function re = RelativeEntropy (p1, p2)

    re = 0;
    for i=1:length(p1)
        if p1(i) > 0.0
            re = re + p1(i)*abs(log(p1(i)/p2(i)));
        end
    end
end


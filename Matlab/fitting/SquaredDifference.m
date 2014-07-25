%  sd = SquaredDifference(p1, p2)
%  
%  Returns the squared difference between two vectors.
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
%  sd : double
%      The squared difference calculated as
%      `sq=\sum_{i=1}^M (p1_i-p2_i)^2`

function sd = SquaredDifference (p1, p2)

    sd = dot(p1-p2, p1-p2);
end


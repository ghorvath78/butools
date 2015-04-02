%  Nm = LagkJointMomentsFromTrace(trace, K, L)
%  
%  Returns the lag-L joint moments of a trace.
%  
%  It is computed as `Nm_{i,j}=\frac{1}{N-L}\sum_{k=0}^{N-L} x_k^i x_{k+L}^j`.
%  
%  Parameters
%  ----------
%  trace : vector of doubles
%      The trace data
%  K : int
%      The joint moments are computed up to order K
%  L : int, optional
%      The lag-L joint moments are computed.
%      The default value is 1.
%  
%  Returns
%  -------
%  Nm : matrix, shape (K,K)
%      The matrix of joint moments, starting from moment 0

function Nm = LagkJointMomentsFromTrace (trace, K, L)
    
    if nargin<3
        L = 1;
    end
    if nargin<2
        K = 3;
    end
    
    Nm = zeros(K,K);
    for i=0:K
        for j=0:K
            Nm(i+1,j+1) = dot(trace(1:end-L).^i, trace(L+1:end).^j) / (length(trace)-L);
        end
    end
end

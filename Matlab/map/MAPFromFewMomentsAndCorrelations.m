%  [D0, D1] = MAPFromFewMomentsAndCorrelations(moms, corr1, r)
%  
%  Creates a Markovian arrival process that has the given
%  2 or 3 marginal moments and lag-1 autocorrelation.
%  The decay of the autocorrelation function can be optionally
%  adjusted as well.
%  The lag-k autocorrelation function `\rho_k` of the 
%  resulting MAP is `\rho_k=r(corr_1/r)^k`.
%  
%  Parameters
%  ----------
%  moms : vector of doubles, length 2 or 3
%      The list of marginal moments to match. 
%  corr1 : double
%      The lag-1 autocorrelation coefficient to match.
%  r : double, optional
%      The decay of the autocorrelation function.
%  
%  Returns
%  -------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the Markovian arrival process
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the Markovian arrival process
%  
%  Notes
%  -----
%  With 2 marginal moments, or with 3 marginal moments and 
%  positive autocorrelation the procedure always returns a 
%  valid Markovian representation.
%  
%  References
%  ----------
%  .. [1] G Horvath, "Matching marginal moments and lag 
%         autocorrelations with MAPs," ValueTools 2013, 
%         Torino, Italy (2013).

function [D0, D1] = MAPFromFewMomentsAndCorrelations (moms, corr1, r)

    m1 = moms(1);
    c2 = moms(2)/moms(1)^2 - 1;
    if length(moms)>2
        l3 = moms(3)*moms(1)/moms(2)^2 - 1;
    else
        l3 = [];
    end

    if ~exist('r','var')
        r = [];
    end
    
    if corr1>=0
        if nargin==5 && ~isempty(r) && r<=0 && r>=1
            error "Parameter r is out of range!";
        end
        if nargin<5 || isempty(r)
            r = 2.0*corr1 / (1.0+corr1);
            p1 = (1.0-(1.0+corr1)/2.0) / (1.0+c2);
            p2 = c2 * p1;
        else
            p1 = (1.0-corr1/r) / (1.0+c2);
            p2 = c2 * p1;
        end
        m11 = m1 * (1.0 - sqrt(r));
        m12 = m1 * (1.0 + c2*sqrt(r));
        if isempty(l3)
            cv21 = (sqrt(c2)*(1.0+c2)*(1.0+sqrt(r))) / (1.0+sqrt(c2)*(1-sqrt(r))+c2*sqrt(r));
            cv22 = (c2*(1.0+c2)*(1-r)) / ((1+c2*sqrt(r)) * (1.0+sqrt(c2)*(1-sqrt(r))+c2*sqrt(r)));
            m21 = (cv21+1.0) * m11^2;
            m22 = (cv22+1.0) * m12^2;
            [alpha1, A1] = APHFrom2Moments ([m11, m21]);
            [alpha2, A2] = APHFrom2Moments ([m12, m22]);
        else
            cv21 = (c2 + sqrt(r)) / (1.0 - sqrt(r));
            cv22 = c2 * (1.0 - sqrt(r)) / (1.0 + c2*sqrt(r));
            l31 = l3 * (c2+1.0) / (c2*(1.0-sqrt(r))+sqrt(c2*(1.0+c2*sqrt(r))*(1.0-sqrt(r))));
            l32 = l3 * (c2+1.0) / ((1.0+c2*sqrt(r))+sqrt(c2*(1.0+c2*sqrt(r))*(1.0-sqrt(r))));        
            m21 = (cv21+1.0) * m11^2;
            m22 = (cv22+1.0) * m12^2;
            m31 = (l31+1.0) * m21^2 / m11;
            m32 = (l32+1.0) * m22^2 / m12;
            [alpha1, A1] = APHFrom3Moments ([m11, m21, m31]);
            [alpha2, A2] = APHFrom3Moments ([m12, m22, m32]);
        end
    else
        if c2>=1
            if nargin==5 && ~isempty(r) && r<=0 && r>=1/c2
                error "Parameter r is out of range!";
            end
            if nargin<5 || isempty(r)
                r = -2.0*corr1 / (1.0-c2*corr1);
                p1 = (1.0+(1.0-c2*corr1)/2.0) / 2.0;
                p2 = p1;
            else
                p1 = (1.0-corr1/r) / 2.0;
                p2 = p1;
            end
        else
            if nargin==5 && ~isempty(r) && r<=0 && r>=1
                error "Parameter r is out of range!";
            end
            if nargin<5 || isempty(r)
                r = -2.0*corr1 / (1.0-corr1);
                p1 = (1.0+(1.0-corr1)/2.0) / 2.0;
                p2 = p1;
            else
                p1 = (1.0-corr1/r) / 2.0;
                p2 = p1;
            end
        end
        m11 = m1 * (1.0 - sqrt(c2*r));
        m12 = m1 * (1.0 + sqrt(c2*r));
        if isempty(l3)
            cv21 = c2 * (1.0-r) / (1.0-sqrt(c2*r));
            cv22 = c2 * (1.0-r) / (1.0+sqrt(c2*r));
            m21 = (cv21+1.0) * m11^2;
            m22 = (cv22+1.0) * m12^2;
            [alpha1, A1] = APHFrom2Moments ([m11, m21]);
            [alpha2, A2] = APHFrom2Moments ([m12, m22]);
        else
            cv21 = (c2+sqrt(c2*r)) / (1.0-sqrt(r*c2));
            cv22 = (c2-sqrt(c2*r)) / (1.0+sqrt(r*c2));
            l31 = 2.0*l3 / (1-sqrt(c2*r)+sqrt(1.0-c2*r));
            l32 = 2.0*l3 / (1+sqrt(c2*r)+sqrt(1.0-c2*r));             
            m21 = (cv21+1.0) * m11^2;
            m22 = (cv22+1.0) * m12^2;
            m31 = (l31+1.0) * m21^2 / m11;
            m32 = (l32+1.0) * m22^2 / m12;
            [alpha1, A1] = APHFrom3Moments ([m11, m21, m31]);
            [alpha2, A2] = APHFrom3Moments ([m12, m22, m32]);
        end
    end

    N1 = length(alpha1);
    N2 = length(alpha2);
    D0 = zeros(N1+N2);
    D0(1:N1,1:N1) = A1;
    D0(N1+1:end,N1+1:end) = A2;
    D1 = zeros(N1+N2);
    D1(1:N1,1:N1) = sum(-A1,2)*alpha1*(1.0-p1);
    D1(1:N1,N1+1:end) = sum(-A1,2)*alpha2*p1;
    D1(N1+1:end,1:N1) = sum(-A2,2)*alpha1*p2;
    D1(N1+1:end,N1+1:end) = sum(-A2,2)*alpha2*(1.0-p2);
end

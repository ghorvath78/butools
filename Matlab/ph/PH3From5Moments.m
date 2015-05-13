%  [alpha, A] = PH3From5Moments(moms)
%  
%  Returns a PH(3) which has the same 5 moments as given.
%  
%  Parameters
%  ----------
%  moms : vector of doubles, length(5)
%    The moments to match
%  
%  Returns
%  -------
%  alpha : vector, shape (1,3)
%      The initial probability vector of the PH(3)
%  A : matrix, shape (3,3)
%      Transient generator matrix of the PH(3)
%  
%  Notes
%  -----
%  Raises an error if the moments are not feasible with
%  a PH(3). Also note that the numerical behavior of the 
%  procedure can be poor if the moments are close to the 
%  boundary of the feasible region.
%  
%  References
%  ----------
%  .. [1] G. Horvath and M. Telek, "On the canonical 
%         representation of phase type distributions," 
%         Performance Evaluation, vol. 66, no. 8, pp. 
%         396 - 409, 2009.

function [alpha, A] = PH3From5Moments (moms, prec)

    if ~exist('prec','var')
        prec = 1e-10;
    end

    % convert the moments to reduced moments
    rmoms = ReducedMomsFromMoms(moms);
    for i=1:5
      rmoms(i) = rmoms(i) / moms(1)^i;
    end

    %solve linear system of equations for a0 a1 a2
    M = [rmoms(3) -rmoms(2) rmoms(1); rmoms(4) -rmoms(3) rmoms(2); rmoms(5) -rmoms(4) rmoms(3)];
    a = M\[1; rmoms(1); rmoms(2)];

    discr = a(3)^2-3*a(2);
    if discr<0
        error ('Invalid characteristic polynomial!');
    end

    gu = (a(3) + 2*sqrt(discr)) / 3;
    g0 = (a(3) + sqrt(discr)) / 3;

    [~,ix] = sort(real(roots([1 fliplr(a')])),'ascend');

    lmb = -roots([1 fliplr(a')]);
    lambda = lmb(ix);

    d1 = a(2) - a(3) - a(1) * rmoms(2);
    d2 = a(1) - a(2) - a(3) * d1;
    d3 = -a(1) - a(2)*d1 - a(3)*d2;

    if d1>prec || (abs(d1)<prec && d2>0)
        error ('Negative density around 0!');
    end

    if lambda(3)<0
        error('Invalid eigenvalues!');
    end

    if imag(lambda(1)) < prec
        gl = real(lambda(1));
    else
        gl = g0;
    end
    
    if gl>gu+prec
        error('Invalid eigenvalues (gl>gu detected)!');
    end
    if gl>gu
        gl = gu;
    end
    
    if abs(d1)<prec
        g2 = 0;
    else
        g2 = -d2 / d1;
    end

    if g2>gu+prec
        error('alpha_2 is negative!');
    end
    if g2>gu
        g2 = gu;
    end

    x1 = max(g2, gl);

    if real(lambda(1))==lambda(1) && g2<gl
        x13 = 0;
    else
        x13 = x1 - a(1) / (x1^2 - a(3)*x1 + a(2));
    end

    bels = (a(3)-x1)^2 - 4*(x1^2-a(3)*x1+a(2));
    if bels<0 && bels>-prec
        bels = 0;
    end

    x2 = (a(3) - x1 + sqrt(bels)) / 2;
    x3 = (a(3) - x1 - sqrt(bels)) / 2;
    p1 = d1 / (x13 - x1);
    p2 = (x1*d1 + d2) / (x13-x1) / x2;
    p3 = (x1*x2*d1 + x2*d2 + x1*d2 + d3) / (x13-x1) / x2 / x3;

    A = [-x1 0 x13; x2 -x2 0; 0 x3 -x3] / moms(1);
    alpha = [p1 p2 p3];

    if x13<-prec || x13>x1
        error('Invalid generator!');
    end

    if min(real(alpha))<-prec
        error('Initial vector has negative entries!');
    end
    
    if max(abs(imag(alpha)))>prec
        error('Inital vector has complex entries!');
    end

    if max(real(alpha))>1+prec
        error('Initial vector has entries that are greater than 1!');
    end
end

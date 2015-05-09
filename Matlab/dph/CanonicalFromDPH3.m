%  [beta, B] = CanonicalFromDPH3(alpha, A, prec)
%  
%  Returns the canonical form of an order-3 discrete phase-type 
%  distribution.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,3)
%      Initial vector of the discrete phase-type distribution
%  A : matrix, shape (3,3)
%      Transition probability matrix of the discrete phase-type
%      distribution
%  prec : double, optional
%    Numerical precision for checking the input, default value
%    is 1e-14
%  
%  Returns
%  -------
%  beta : matrix, shape (1,3)
%    The initial probability vector of the canonical form
%  B : matrix, shape (3,3)
%    Transition probability matrix of the canonical form

function [beta,B] = CanonicalFromDPH3 (alpha,A,prec)
   
    if ~exist('prec','var')
        prec = 1e-14;
    end

    global BuToolsCheckInput;

    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMGRepresentation(alpha,A)
        error('CanonicalFromDPH3: Input isn''t a valid MG distribution!');
    end

    if BuToolsCheckInput && (size(A,1)~=3 || size(A,2)~=3)
        error('CanonicalFromDPH3: Dimension must be 3!');
    end

    lambda = EigSort(eig(A));
    a0 = -lambda(1)*lambda(2)*lambda(3);
    a1 = lambda(1)*lambda(2)+lambda(1)*lambda(3)+lambda(2)*lambda(3);
    a2 = -lambda(1)-lambda(2)-lambda(3);
    e=[1;1;1];

    if real(lambda(1))>0 && real(lambda(2))>=0 && real(lambda(3))>=0
        %PPP case
        [alphaout,A2] = CanonicalFromPH3(alpha,A-eye(3),prec);
        Aout = A2+eye(3);
    elseif real(lambda(1))>0 && real(lambda(2))>=0 && real(lambda(3))<0
        %PPN case
        x1 = lambda(1);
        x2 = lambda(2)+lambda(3);
        x3=lambda(2)*lambda(3)/(lambda(2)+lambda(3)-1);
        Aout=[x1,1-x1,0;0,x2,1-x2;0,x3,0];
        b3=1/(1-x3)*(e-A*e);
        b2=1/(1-x2)*A*b3;
        b1=e-b2-b3;
        B=[b1,b2,b3];
        alphaout=alpha*B;
    elseif real(lambda(1))>0 && real(lambda(2))<0 && real(lambda(3))>=0
        %PNP case
        x1=-a2;
        x2=(a0-a1*a2)/(a2*(1+a2));
        x3=a0*(1+a2)/(a0-a2-a1*a2-a2^2);
        Aout=[x1,1-x1,0;x2,0,1-x2;0,x3,0];
        b3=1/(1-x3)*(e-A*e);
        b2=1/(1-x2)*A*b3;
        b1=e-b2-b3;
        if alpha*b1>=0
            B=[b1,b2,b3];
            alphaout = alpha*B;
        else
            %Set the initial vector first element to 0
            x33=0;
            a1=-1;
            while x33 <= 1
                [a1,x1,x2,x3,B]=firstInitElem(x33,lambda,alpha,A);
                if a1 >= 0 && x1 >= 0 && x2 >= 0 && x3 >= 0 && x3+x33 < 1
                    break;
                end
                x33=x33+0.01;
            end
            
            if a1 >= 0
                Aout=[x1,1-x1,0;x2,0,1-x2;0,x3,x33];
                alphaout=alpha*B;
            else
                %PNP+
                x1=lambda(3);
                x2=lambda(1)+lambda(2);
                x3=lambda(1)*lambda(2)/(lambda(1)+lambda(2)-1);
                Aout=[x1,0,0;0,x2,1-x2;0,x3,0];

                p1=alpha*(e-A*e);
                p2=alpha*A*(e-A*e);
                d1=(1-lambda(1))*((1-lambda(2))*(1-lambda(3))+(-1+lambda(2)+lambda(3))*p1-p2)/((lambda(1)-lambda(2))*(lambda(1)-lambda(3)));
                d2=(lambda(2)-1)*((1-lambda(1))*(1-lambda(3))+(-1+lambda(1)+lambda(3))*p1-p2)/((lambda(1)-lambda(2))*(lambda(2)-lambda(3)));
                d3=(lambda(3)-1)*((1-lambda(1))*(1-lambda(2))+(-1+lambda(1)+lambda(2))*p1-p2)/((lambda(2)-lambda(3))*(lambda(3)-lambda(1)));
                alphaout={d3/(1-lambda(3)),(d1*lambda(1)+d2*lambda(2))/((1-lambda(1))*(1-lambda(2))),(d1+d2)*(1-lambda(1)-lambda(2))/((1-lambda(1))*(1-lambda(2)))}; 

			    if min(alphaout) < 0 || min(min(Aout)) < 0
				    fprintf('DPH3Canonical: Unhandled PNP case! Input:\n');
				    disp(alpha);
				    disp(A);
				    error('DPH3Canonical: PNP ERROR!');
			    end
            end
        end
    elseif real(lambda(1))>0 && real(lambda(2))<0 && real(lambda(3))<0
        %PNN case
        if isreal(lambda) || abs(lambda(2))^2 <= 2*lambda(1)*(-real(lambda(2)))
            x1=-a2;
            x2=-a1/(1+a2);
            x3=-a0/(1+a1+a2);
            Aout=[x1,1-x1,0;x2,0,1-x2;x3,0,0];
            b3=1/(1-x3)*(e-A*e);
            b2=1/(1-x2)*A*b3;
            b1=e-b2-b3;
            B=[b1,b2,b3];
            alphaout=alpha*B;
        else
           [alphaout,A2]=PH3Canonical(alpha,A-eye(3));
           Aout=A2+eye(3);
        end
    end
    beta = alphaout;
    B = Aout;
end

function [a1,m1,m2,m3,B]=firstInitElem(m33,sortEigs,alpha,A)
    l1=sortEigs(1);
    l2=sortEigs(2);
    l3=sortEigs(3);
    
    m1=-m33+l1+l2+l3;
    m2=-((l2-l3)*(l1^2-l1*l2-l1*l3+l2*l3)*(m33^3-2*m33^2*l1+m33*l1^2-2*m33^2*l2+3*m33*l1*l2-l1^2*l2+m33*l2^2- ...
      l1*l2^2-2*m33^2*l3+3*m33*l1*l3-l1^2*l3+3*m33*l2*l3-2*l1*l2*l3-l2^2*l3+m33*l3^2-l1*l3^2-l2*l3^2))/ ...
      (2*m33*l1^2*l2+2*m33^2*l1^2*l2-l1^3*l2-3*m33*l1^3*l2+l1^4*l2-2*m33*l1*l2^2-2*m33^2*l1*l2^2+l1^3*l2^2+ ...
      l1*l2^3+3*m33*l1*l2^3-l1^2*l2^3-l1*l2^4-2*m33*l1^2*l3-2*m33^2*l1^2*l3+l1^3*l3+3*m33*l1^3*l3-l1^4*l3+ ...
      2*m33*l2^2*l3+2*m33^2*l2^2*l3-l2^3*l3-3*m33*l2^3*l3+l2^4*l3+2*m33*l1*l3^2+2*m33^2*l1*l3^2-l1^3*l3^2- ...
      2*m33*l2*l3^2-2*m33^2*l2*l3^2+l2^3*l3^2-l1*l3^3-3*m33*l1*l3^3+l1^2*l3^3+l2*l3^3+3*m33*l2*l3^3-l2^2*l3^3+ ...
      l1*l3^4-l2*l3^4);
    m3=((l2-l3)*(l1^2-l1*l2-l1*l3+l2*l3)*(m33^3+m33^4-m33^2*l1-2*m33^3*l1+m33^2*l1^2-m33^2*l2-2*m33^3*l2+ ...
      m33*l1*l2+3*m33^2*l1*l2-m33*l1^2*l2+m33^2*l2^2-m33*l1*l2^2-m33^2*l3-2*m33^3*l3+m33*l1*l3+3*m33^2*l1*l3- ...
      m33*l1^2*l3+m33*l2*l3+3*m33^2*l2*l3-l1*l2*l3-4*m33*l1*l2*l3+l1^2*l2*l3-m33*l2^2*l3+l1*l2^2*l3+m33^2*l3^2- ...
      m33*l1*l3^2-m33*l2*l3^2+l1*l2*l3^2))/(-2*m33*l1^2*l2-2*m33^2*l1^2*l2-m33^3*l1^2*l2+l1^3*l2+3*m33*l1^3*l2+ ...
      2*m33^2*l1^3*l2-l1^4*l2-m33*l1^4*l2+2*m33*l1*l2^2+2*m33^2*l1*l2^2+m33^3*l1*l2^2-l1^3*l2^2-2*m33*l1^3*l2^2+ ...
      l1^4*l2^2-l1*l2^3-3*m33*l1*l2^3-2*m33^2*l1*l2^3+l1^2*l2^3+2*m33*l1^2*l2^3+l1*l2^4+m33*l1*l2^4-l1^2*l2^4+ ...
      2*m33*l1^2*l3+2*m33^2*l1^2*l3+m33^3*l1^2*l3-l1^3*l3-3*m33*l1^3*l3-2*m33^2*l1^3*l3+l1^4*l3+m33*l1^4*l3- ...
      2*m33*l2^2*l3-2*m33^2*l2^2*l3-m33^3*l2^2*l3+l2^3*l3+3*m33*l2^3*l3+2*m33^2*l2^3*l3-l2^4*l3-m33*l2^4*l3- ...
      2*m33*l1*l3^2-2*m33^2*l1*l3^2-m33^3*l1*l3^2+l1^3*l3^2+2*m33*l1^3*l3^2-l1^4*l3^2+2*m33*l2*l3^2+2*m33^2*l2*l3^2+ ...
      m33^3*l2*l3^2-l2^3*l3^2-2*m33*l2^3*l3^2+l2^4*l3^2+l1*l3^3+3*m33*l1*l3^3+2*m33^2*l1*l3^3-l1^2*l3^3- ...
      2*m33*l1^2*l3^3-l2*l3^3-3*m33*l2*l3^3-2*m33^2*l2*l3^3+l2^2*l3^3+2*m33*l2^2*l3^3-l1*l3^4-m33*l1*l3^4+ ...
      l1^2*l3^4+l2*l3^4+m33*l2*l3^4-l2^2*l3^4);
    b3=sum(eye(3)-A,2)/(1-m3-m33);
    b2=(-m33*eye(3)+A)*b3/(1-m2);
    b1=(-m3*b3+A*b2)/(1-m1);
    B=[b1,b2,b3];
    a1=alpha*b1;
end

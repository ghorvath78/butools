%  jm = JMomsFromJFactorialMoms(jfm)
%  
%  Returns the lag-1 joint raw moments given the lag-1 joint
%  factorial moments.
%  
%  The lag-1 joint raw moments are: 
%  `m_{i,j}=E(\mathcal{X}^i\mathcal{Y}^j)`
%  
%  The factorial moments are: 
%  `f_{ij}=E(\mathcal{X}(\mathcal{X}-1)\cdots(\mathcal{X}-i+1)\mathcal{Y}(\mathcal{Y}-1)\cdots(\mathcal{Y}-j+1))`
%     
%  Parameters
%  ----------
%  jfm : matrix, shape (M,M)
%      The matrix of joint factorial moments. The entry in 
%      row i and column j is `f_{i,j},i\geq 1,j\geq 1`.
%      
%  Returns
%  -------
%  jm : matrix, shape (M,M)
%      The matrix of joint raw moments. The entry in row i
%      and column j is `m_{i,j},i\geq 1,j\geq 1`.
%      
%  References
%  ----------
%  http://en.wikipedia.org/wiki/Factorial_moment    

function jmoms=JMomsFromJFactorialMoms(jfmoms)

    s1=size(jfmoms,1);
    s2=size(jfmoms,2);

    jmoms=zeros(s1,s2);
    for i=1:s1
        for j=1:s2
            xCoeff=poly(0:i-1);
            xCoeff=xCoeff(end-1:-1:1);
            yCoeff=poly(0:j-1);
            yCoeff=yCoeff(end-1:-1:1);
            eh=-xCoeff'*yCoeff;
            jmoms(i,j)=jfmoms(i,j)+trace(jmoms(1:i,1:j)*eh');
        end
    end
end
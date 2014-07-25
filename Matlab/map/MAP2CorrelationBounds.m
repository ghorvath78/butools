%  [lb, ub] = MAP2CorrelationBounds(moms)
%  
%  Returns the upper and lower correlation bounds for a MAP(2)
%  given the three marginal moments.
%  
%  !!!TO CHECK!!!
%  
%  Parameters
%  ----------
%  moms : vector, length(3)
%      First three marginal moments of the inter-arrival times
%  
%  Returns
%  -------
%  lb : double
%      Lower correlation bound
%  ub : double
%      Upper correlation bound
%  
%  References
%  ----------
%  .. [1] L Bodrog, A Heindl, G Horvath, M Telek, "A Markovian
%         Canonical Form of Second-Order Matrix-Exponential 
%         Processes," EUROPEAN JOURNAL OF OPERATIONAL RESEARCH
%         190:(2) pp. 459-477. (2008)
%         

function [lb, ub] = MAP2CorrelationBounds (moms)

    m1 = moms(1);
    m2 = moms(2);
    m3 = moms(3);

    h2 = m2/(2.0*m1*m1) - 1;
    h3 = m3/(6.0*m1*m1*m1)-m2*m2/(4.0*m1*m1*m1*m1);
    cv2 = m2/m1/m1 - 1.0;

    if h2>=0  
        gub = h2;
    else
        gub = -(h2+sqrt(-h3))^2;
    end

    if h2<=0 || h3/h2+h2<1
        glb = -h3 - h2*h2;
    else
        glb = h2 * (h3+h2*h2-h2-sqrt((h3+h2*h2-h2)^2+4.0*h2*h2*h2)) / (h3+h2*h2-h2+sqrt((h3+h2*h2-h2)^2+4.0*h2*h2*h2));
    end
    
    if h2>=0
        lb = glb/cv2;
        ub = gub/cv2;
    else
        ub = glb/cv2;
        lb = gub/cv2;
    end
end

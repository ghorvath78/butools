%  D = RandomDMMAP(order, types, mean, zeroEntries, maxTrials, prec)
%  
%  Returns a random discrete Markovian arrival process.
%  
%  Parameters
%  ----------
%  order : int
%      The size of the DMAP
%  mean : double, optional
%      The mean inter-arrival times of the DMMAP
%  types : int
%      The number of different arrival types
%  zeroEntries : int, optional
%      The number of zero entries in the D0 and D1 matrices
%  maxTrials : int, optional
%      The maximum number of trials to find a proper DMMAP 
%      (that has an irreducible phase process and none of 
%      its parameters is all-zero)
%  prec : double, optional
%      Numerical precision for checking the irreducibility.
%      The default value is 1e-14.
%  
%  Returns
%  -------
%  D : list/cell of matrices of shape(M,M), length(types+1)
%      The D0...Dtypes matrices of the DMMAP 
%  
%  Notes
%  -----
%  If it fails, try to increase the 'maxTrials' parameter,
%  or/and the 'mean' parameter.

function D = RandomDMMAP(order, types, mean, zeroEntries, maxTrials, prec)
% RandomDMMAP [ order, types, zeroEntries[0], prec[10^-14], maxTrials[100] ] 
%  -> [ cell of matrix0, matrix1 .. matrixM ] :
%     Generates a random discrete time marked Markovian arrival process of 
%     the given order and the given number of types. The obtained 
%     representation containes 'zeroEntries' zeros. If it fails after 
%     'maxTrials' trials, then it decreases the number of zero entries. It
%     prints a message, if the found representation contains less zeros 
%     than given. 'prec' is the numerical precision

    if ~exist('zeroEntries','var')
        zeroEntries = 0;
    end

    if ~exist('mean','var')
        mean = 10;
    end

    if ~exist('prec','var')
        prec = 1e-7;
    end

    if ~exist('maxTrials','var')
        maxTrials = 1000;
    end

    if types < 1
        error('RandomMMAP: ''types'' must be positive integer!');
    end

    if zeroEntries > (order+1)*(order-1)+types*(order^2-1)
        error('RandomDMAP/DMMAP: Too many zero entries requested! Try to decrease the zeroEntries parameter!');
    end

    % distribute the zero entries among the rows
    function o = allZeroDistr (states, zeros)
        if states==1
            o = zeros;
        else
            o = [];
            for iz=0:zeros
                x = allZeroDistr (states-1, zeros-iz);
                for jz=1:size(x,1)
                    xt = sort([x(jz,:), iz]);
                    % check if we have it already
                    found = 0;
                    for kz=1:size(o,1)
                        if sum((o(kz,:)-xt).^2)==0
                            found = 1;
                            break;
                        end
                    end
                    if ~found
                        o = [o; xt];
                    end
                end
            end
        end
    end
        
    zeroDistr = allZeroDistr(order, zeroEntries);   

    trials = 1;
    while trials<maxTrials
        % select a configuration from zeroDistr: it is a list describing the zero entries in each row
        zdix = randperm(size(zeroDistr,1));
        for k=1:size(zeroDistr,1)
            zDistr = zeroDistr(zdix(k),:); 
            bad = 0;
            for di=1:length(zDistr)
                if zDistr(di)>=(types+1)*order-1;
                    bad = 1;
                    break;
                end
            end
            if bad
                trials = trials + 1;
                continue;
            end            
            B = zeros(order,(types+1)*order);
            for i=1:order
                rp = randperm((types+1)*order-1);
                a = zeros(1,(types+1)*order-1);
                for j=1:(types+1)*order-1-zDistr(i)
                    a(rp(j)) = rand();
                end
                B(i,1:i-1) = a(1:i-1);
                B(i,i+1:end) = a(i:end);
            end           
            % construct DMMAP matrices
            D = cell(1,types+1);
            sc = zeros(order,1);
            for i=1:types+1
                D{i} = B(:,(i-1)*order+1:i*order);
                sc = sc + sum(D{i},2);
            end
            if any(sc==0)
                continue;
            end
            for i=1:types+1
                D{i} = diag(1./sc)*D{i};
            end
            
            % check if it is a proper MAP (irreducible phase process & no full zero matrix)
            sumD = D{1};
            for i=2:length(D)
                sumD = sumD + D{i};
            end
            if rank(D{1})==order && rank(eye(order)-sumD)==order-1
                alpha = DTMCSolve(sumD);
                if min(abs(alpha)) > prec
                    fullZero = 0;
                    for i=1:length(D)
                        if all(all(D{i}==0.0))
                            fullZero = 1;
                            break;
                        end
                    end
                    if fullZero==0
                        % diagonals of matrix D0:
                        d = rand(1,order);
                        % scale to the mean value
                        Dv = cell(1,types+1);
                        for i=1:types+1
                            Dv{i} = diag(1-d)*D{i};
                        end
                        Dv{1} = Dv{1}+diag(d);                        
                        m = MarginalMomentsFromDMMAP (Dv, 1);
                        d = 1 - (1-d)*m(1)/mean;
                        for i=1:types+1
                            D{i} = diag(1-d)*D{i};
                        end
                        D{1} = D{1}+diag(d);                        
                        if CheckDMMAPRepresentation(D)
                            return;
                        end
                    end
                end
            end
            trials = trials + 1;
        end
    end
    error('No feasible random DMAP/DMMAP found with such many zero entries! Try to increase the maxTrials parameter, or the mean value!');
end


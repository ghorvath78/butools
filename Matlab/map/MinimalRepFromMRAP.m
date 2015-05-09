%  D = MinimalRepFromMRAP(H, how, precision)
%  
%  Returns the minimal representation of a marked rational
%  arrival process.
%  
%  Parameters
%  ----------
%  H : list of matrices of shape (M,M)
%      The list of H0, H1, ..., HK matrices of the marked
%      rational arrival process
%  how : {"obs", "cont", "obscont"}, optional        
%      Determines how the representation is minimized. 
%      "cont" means controllability, "obs" means 
%      observability, "obscont" means that the rational arrival
%      process is minimized in both respects. Default value 
%      is "obscont".
%  precision : double, optional
%     Precision used by the Staircase algorithm. The default
%     value is 1e-12.
%  
%  Returns
%  -------
%  D : list of matrices of shape (M,M)
%      The D0, D1, ..., DK matrices of the minimal 
%      representation
%  
%  References
%  ----------
%  .. [1] P. Buchholz, M. Telek, "On minimal representation of 
%         rational arrival processes." Madrid Conference on 
%         Qeueuing theory (MCQT), June 2010.

function D = MinimalRepFromMRAP (H, how, precision)

    if ~exist('precision','var')
        precision = 1e-12;
    end

    if ~exist('how','var')
        how = 'obscont';
    end
    
    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end       

    if BuToolsCheckInput && ~CheckMRAPRepresentation (H)
        error('MinimalRepFromMRAP: Input isn''t a valid MRAP representation');
    end

    if strcmp(how,'cont')
        [B, n] = MStaircase (H, ones(size(H{1},1),1), precision);
        D = cell(1,length(H));
        for i=1:length(H)
            Di = inv(B)*H{i}*B;
            D{i} = Di(1:n,1:n);
        end
    elseif strcmp(how,'obs')
        [alpha, A] = MarginalDistributionFromMRAP (H);
        G = cell(1,length(H));
        for i=1:length(H)
            G{i} = H{i}';
        end
        [B, n] = MStaircase (G, alpha', precision);
        for i=1:length(H)
            Di = inv(B)*H{i}*B;
            D{i} = Di(1:n,1:n);
        end
    elseif strcmp(how,'obscont')
        D = MinimalRepFromMRAP(H, 'cont', precision);
        D = MinimalRepFromMRAP(D, 'obs', precision);
    end
end


%  [D0, D1] = MinimalRepFromRAP(H0, H1, how, precision)
%  
%  Returns the minimal representation of a rational arrival 
%  process.
%  
%  Parameters
%  ----------
%  H0 : matrix, shape (M,M)
%      The H0 matrix of the rational arrival process
%  H1 : matrix, shape (M,M)
%      The H1 matrix of the rational arrival process
%  how : {"obs", "cont", "obscont"}, optional      
%      Determines how the representation is minimized. "cont" 
%      means controllability, "obs" means observability, 
%      "obscont" means that the rational arrival process is 
%      minimized in both respects. The default value is 
%      "obscont"
%  precision : double, optional
%     Precision used by the Staircase algorithm. The default 
%     value is 1e-12.
%  
%  Returns
%  -------
%  D0 : matrix, shape (M,M)
%      The D0 matrix of the minimal representation
%  D1 : matrix, shape (M,M)
%      The D1 matrix of the minimal representation
%  
%  References
%  ----------
%  .. [1] P. Buchholz, M. Telek, "On minimal representation of 
%         rational arrival processes." Madrid Conference on 
%         Qeueuing theory (MCQT), June 2010.

function [D0, D1] = MinimalRepFromRAP (H0, H1, how, precision)

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

    if BuToolsCheckInput && ~CheckRAPRepresentation (H0, H1)
        error('MinimalRepFromRAP: Input isn''t a valid MRAP representation');
    end

    D = MinimalRepFromMRAP ({H0,H1}, how, precision);
    D0 = D{1};
    D1 = D{2};
end


%  B = TransformToAcyclic(A, maxSize, precision)
%  
%  Transforms an arbitrary matrix to a Markovian bi-diagonal 
%  matrix.
%  
%  Parameters
%  ----------
%  A : matrix, shape (N,N)
%      Matrix parameter of the initial representation
%  maxSize : int, optional
%      The maximal order of the resulting Markovian 
%      representation. The default value is 100
%  precision : double, optional
%      Matrix entries smaller than the precision are 
%      considered to be zeros. The default value is 1e-14
%      
%  Returns
%  -------
%  B : matrix, shape (N,N)
%      Transient (bi-diagonal) generator matrix of the
%      Markovian acyclic representation.
%      
%  Notes
%  -----
%  Calls the 'TransformToMonocyclic' procedure if all the 
%  eigenvalues are real, otherwise it raises an error if no
%  Markovian acyclic generator has been found up to the 
%  given size.
%  
%  Raises an error if no Markovian acyclic generator 
%  has been found up to the given size.
%  %  References
%  ----------
%  .. [1]  Mocanu, S., Commault, C.: "Sparse representations 
%          of phase-type distributions," Stoch. Models 15, 
%          759-778 (1999)

function B = TransformToAcyclic (A, maxSize, precision)

    if ~exist('precision','var')
        precision = 1e-14;
    end

    if ~exist('maxSize','var')
        maxSize = 100;
    end

    if any(imag(eig(A))>=precision)
        error 'TransformToAcyclic: Complex eigenvalue found, no acyclic representation exists.';
    end

    B = TransformToMonocyclic (A, maxSize, precision);
end

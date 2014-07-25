%  ImageFromDPH(alpha, A, outFileName, prec)
%  
%  Depicts the given discrete phase-type distribution,
%  and either displays it or saves it to file.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,M)
%      The initial probability vector of the discrete phase-
%      type distribution.
%  A : matrix, shape (M,M)
%      The transition probability  matrix of the discrete phase-
%      type distribution.
%  outFileName : string, optional
%      If it is not provided, or equals to 'display', the
%      image is displayed on the screen, otherwise it is 
%      written to the file. The file format is deduced 
%      from the file name.
%  prec : double, optional
%      Transition probabilities less then prec are 
%      considered to be zero and are left out from the 
%      image. The default value is 1e-13.
%  
%  Notes
%  -----
%  The 'graphviz' software must be installed and available
%  in the path to use this feature.

function img = ImageFromPH(alpha,A,outFileName,prec)

    if ~exist('prec','var')
        prec = 1e-13;
    end

    if ~exist('outFileName','var') || strcmp(outFileName,'display')
        outputFile = '.result.png';
        displ = true;
    else
        outputFile = outFileName;
        displ = false;
    end
    
    inputFile = '.temp.dot';

    fid = fopen(inputFile,'w');
    fprintf(fid, 'digraph G {\n');
    fprintf(fid, '\trankdir=LR;\n');
    fprintf(fid, '\tnode [shape=circle,width=0.3,height=0.3,label=""];\n');
    
    % nodes
    for i=1:length(alpha)
        fprintf(fid, '\tn%d [xlabel=<<i>%g</i>>];\n', i, alpha(i));
    end
    
    % transitions to a non-absorbing state
    for i=1:length(alpha)
        for j=1:length(alpha)
            if abs(A(i,j))>prec
                fprintf(fid, '\tn%d -> n%d [label="%g"];\n', i, j, A(i,j));
            end
        end
    end
    
    % transitions to the absorbing state
    fprintf(fid, ['\tab [style=filled];\n']);    
    a = 1-sum(A,2);
    for i=1:length(alpha)
        if abs(a(i))>prec
            fprintf(fid, '\tn%d -> ab [label="%g"];\n', i, a(i));
        end
    end    
    fprintf(fid,'}\n');
    fclose(fid);

    [~,~,ext] = fileparts(outputFile);
    system(['dot -T', ext(2:end), ' ', inputFile, ' -o ', outputFile]);
    
    delete (inputFile);

    if displ
        RGB = imread(outputFile);
        figure('toolbar','none','units','pixel','position',[0,0,size(RGB,2)+100,size(RGB,1)+100]);
        image(RGB);  
        axis image;
        delete(outputFile);
    end
end


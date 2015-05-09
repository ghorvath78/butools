%  ImageFromMMAP(D, outFileName, prec)
%  
%  Depicts the given marked Markovian arrival process, and
%  either displays it or saves it to file.
%  
%  Parameters
%  ----------
%  D : list of matrices of shape(M,M), length(N)
%      The D0...DN matrices of the MMAP
%  outFileName : string, optional
%      If it is not provided, or equals to 'display', the
%      image is displayed on the screen, otherwise it is 
%      written to the file. The file format is deduced 
%      from the file name.
%  prec : double, optional
%      Transition rates less then prec are considered to
%      be zero and are left out from the image. The 
%      default value is 1e-13.
%  
%  Notes
%  -----
%  The 'graphviz' software must be installed and available
%  in the path to use this feature.

function ImageFromMMAP(D,outFileName,prec)

    global BuToolsCheckInput;
    if isempty(BuToolsCheckInput)
        BuToolsCheckInput = true;
    end   

    if BuToolsCheckInput && ~CheckMAPRepresentation(D)
        error('ImageFromMMAP: Input isn''t a valid MMAP representation!');
    end
    
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

    N = size(D{1},1);
    
    % transitions without arrivals
    Dx=D{1};
    for i=1:N
        for j=1:N
            if i~=j && abs(Dx(i,j))>prec
                fprintf(fid, '\tn%d -> n%d [label="%g"];\n', i, j, Dx(i,j));
            end
        end
    end
    
    % transitions with arrivals
    for k=2:length(D)
        Dx=D{k};
        for i=1:N
            for j=1:N
                if abs(Dx(i,j))>prec
                    if length(D)==2
                        fprintf(fid, '\tn%d -> n%d [style="dashed",label="%g"];\n', i, j, Dx(i,j));
                    else
                        fprintf(fid, '\tn%d -> n%d [style="solid",fontcolor="/dark28/%d",color="/dark28/%d",label="%g"];\n', i, j, min(k-1,8), min(k-1,8), Dx(i,j));
                    end
                end
            end
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


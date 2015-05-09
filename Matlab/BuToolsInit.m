function BuToolsInit (verbose, checkInput, checkPrecision)

    fname = mfilename('fullpath');
    [pathstr, ~, ~] = fileparts(fname);

    disp('Butools V2.0')
    
    packages = {'utils', 'mc', 'moments', 'reptrans', 'trace', 'ph', 'dph', 'map', 'dmap', 'mam', 'queues', 'fitting'};
    
    path(pathstr, path);
    fprintf('Packages loaded: ');
    for p=1:length(packages)
        path(fullfile(pathstr,packages{p}), path);
        fprintf(packages{p});
        if p<length(packages)
            fprintf(', ');
        else
            fprintf('\n');
        end
    end   
    
    global BuToolsVerbose;
    global BuToolsCheckInput;
    global BuToolsCheckPrecision;
    
    if exist('verbose','var')
        BuToolsVerbose = verbose;
    else
        BuToolsVerbose = false;
    end

    if exist('checkInput','var')
        BuToolsCheckInput = checkInput;
    else
        BuToolsCheckInput = true;
    end
    
    if exist('checkPrecision','var')
        BuToolsCheckPrecision = checkPrecision;
    else
        BuToolsCheckPrecision = 1e-12;
    end
    
    fprintf('Global variables: BuToolsVerbose = %d, BuToolsCheckInput = %d, BuToolsCheckPrecision = %g\n',BuToolsVerbose,BuToolsCheckInput,BuToolsCheckPrecision);
end


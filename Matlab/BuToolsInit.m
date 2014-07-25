function BuToolsInit (verbose, checkInput)

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
    
    fprintf('Global variables: BuToolsVerbose = %d, BuToolsCheckInput = %d\n',BuToolsVerbose,BuToolsCheckInput);   
end


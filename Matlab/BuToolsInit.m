function BuToolsInit (verbose, checkInput, checkPrecision)

    fname = mfilename('fullpath');
    [pathstr, ~, ~] = fileparts(fname);

    disp('<strong>Butools V2.0</strong>')
    
    packages = {'utils', 'mc', 'moments', 'reptrans', 'trace', 'ph', 'dph', 'map', 'dmap', 'mam', 'queues', 'fitting'};
    
    pkgpath = {pathstr};
    for p=1:length(packages)
        pkgpath{p+1} = fullfile(pathstr,packages{p});
    end    
    addpath(strjoin(pkgpath,':'));
    
    fprintf('Packages loaded: ');
    for p=1:length(packages)
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


%%
% surrogate command line tool

function surrogate(varargin)

    % set version number
    versionNumber = '0.1';

    % add script path
    if ~isdeployed % checking MATLAB mode or stand-alone mode.
        [st,ind] = dbstack('-completenames');
        relpath = st(ind).file;
        [exedir,exename,ext] = fileparts(relpath);
        if exist([exedir '/util'],'dir')
            addpath([exedir '/util']);
            addpath([exedir '/lib']);
        end
    end

    % get exe file full path
    global exePath;
    global exeName;
    [exePath, exeName, ext] = exeFilename();

    % init command line input
    handles.commandError = 0;
    handles.csvFiles = {};
    handles.outpath = 'results';
    handles.lag = 3;
    handles.transform = 0;
    handles.transopt = NaN;
    handles.format = 0;
    handles.showSig = 0;

    handles.rg = 0;
    handles.rs = 0;
    handles.ft = 0;
    handles.aaft = 0;
    handles.iaaft = 0;
    handles.var = 0;
    handles.pcvar = 0;
    handles.vardnn = 0;
    handles.lazy = 0;
    handles.multi = 1;
    handles.uni = 0;
    handles.noiseType = 'gaussian';
    handles.surrNum = 1;

    handles.maxEpochs = 1000;
    handles.L2Regularization = 0.05;
    handles.noCache = 0;
    handles.nn = 2;

    % load command line input
    i = 1;
    while true
        if i > size(varargin, 2)
            break;
        end
        switch varargin{i}
            case {'-g','--rg'}
                handles.rg = 1;
            case {'-s','--rs'}
                handles.rs = 1;
            case {'-f','--ft'}
                handles.ft = 1;
            case {'-a','--aaft'}
                handles.aaft = 1;
            case {'-i','--iaaft'}
                handles.iaaft = 1;
            case {'-v','--var'}
                handles.var = 1;
            case {'-p','--pcvar'}
                handles.pcvar = 1;
            case {'-d','--vardnn'}
                handles.vardnn = 1;
            case {'-l','--lazy'}
                handles.lazy = 1;
            case {'--multi'}
                handles.multi = 1;
            case {'--uni'}
                handles.uni = 1;
            case {'--noise'}
                handles.noiseType = varargin{i+1};
                i = i + 1;
            case {'--outnum'}
                handles.surrNum = str2num(varargin{i+1});
                i = i + 1;
            case {'--outpath'}
                handles.outpath = varargin{i+1};
                i = i + 1;
            case {'--lag'}
                handles.lag = str2num(varargin{i+1});
                i = i + 1;
            case {'--transform'}
                handles.transform = str2num(varargin{i+1});
                i = i + 1;
            case {'--transopt'}
                handles.transopt = str2num(varargin{i+1});
                i = i + 1;
            case {'--format'}
                handles.format = str2num(varargin{i+1});
                i = i + 1;
            case {'--epoch'}
                handles.maxEpochs = str2num(varargin{i+1});
                i = i + 1;
            case {'--l2'}
                handles.L2Regularization = str2num(varargin{i+1});
                i = i + 1;
            case {'--nn'}
                handles.nn = str2num(varargin{i+1});
                i = i + 1;
            case {'--nocache'}
                handles.noCache = 1;
            case {'--showsig'}
                handles.showSig = 1;
            case {'-h','--help'}
                showUsage();
                return;
            case {'-v','--version'}
                disp([exeName ' version : ' num2str(versionNumber)]);
                return;
            otherwise
                if strcmp(varargin{i}(1), '-')
                    disp(['bad option : ' varargin{i}]);
                    i = size(varargin, 2);
                    handles.commandError = 1;
                else
                    handles.csvFiles = [handles.csvFiles varargin{i}];
                end
        end
        i = i + 1;
    end
    
    % check command input
    if handles.commandError
        showUsage();
        return;
    elseif isempty(handles.csvFiles)
        disp('no input files. please specify node status signal files.');
        showUsage();
        return;
    end

    % process input files
    processInputFiles(handles);
end

%%
% show usage function
function showUsage()
    global exePath;
    global exeName;
    
    disp(['usage: ' exeName ' [options] filename.csv ...']);
    disp('  -g, --rg            output Random Gaussian (RG) surrogate (<filename>_rg_<variate>_<num>.csv)');
    disp('  -s, --rs            output Random Shuffling (RS) surrogate (<filename>_rs_<variate>_<num>.csv)');
    disp('  -f, --ft            output Fourier Transform (FT) surrogate (<filename>_ft_<variate>_<num>.csv)');
    disp('  -a, --aaft          output Amplitude Adjusted FT (AAFT) surrogate (<filename>_aaft_<variate>_<num>.csv)');
    disp('  -i, --iaaft         output Iterated AAFT (IAAFT) surrogate (<filename>_iaaft_<variate>_<num>.csv)');
    disp('  -v, --var           output Vector Auto-Regression (VAR) surrogate (<filename>_var_<variate>_<num>.csv)');
    disp('  -p, --pcvar         output Principal Component VAR (PCVAR) surrogate (<filename>_pcvar_<variate>_<num>.csv)');
    disp('  -d, --vardnn        output VAR Deep Neural Network (VARDNN) surrogate (<filename>_vardnn_<variate>_<num>.csv)');
    disp('  -l, --lazy          output Lazy Learning (LL) surrogate (<filename>_lazy_<variate>_<num>.csv)');
    disp('  --multi             output multivariate surrogate (default:on)');
    disp('  --uni               output univariate surrogate (default:off)');
    disp('  --noise type        noise type for VAR, PCVAR, VARDNN, LL surrogate (default:"gaussian")');
    disp('  --outnum num        output surrogate sample number <num> (default:1)');
    disp('  --outpath           output files path (default:"results")');
    disp('  --format type       save file format <type> 0:csv, 1:mat(each), 2:mat(all) (default:0)');
    disp('  --transform type    input signal transform <type> 0:raw, 1:sigmoid (default:0)');
    disp('  --transopt num      signal transform option <num> (for type 1:centroid value)');
    disp('  --lag num           time lag <num> for VAR, PCVAR, VARDNN, LL (default:3)');
    disp('  --epoch num         VARDNN training epoch number <num> (default:1000)');
    disp('  --l2 num            VARDNN training L2Regularization <num> (default:0.05)');
    disp('  --nn num            <num>-nearest neighbor for Lazy Learning (default:2)');
    disp('  --showsig           show node status signals of <filename>.csv');
    disp('  --nocache           do not use cache file for VARDNN training');
    disp('  -v, --version       show version number');
    disp('  -h, --help          show command line help');
end

%%
% process input files (mail rutine)
%
function processInputFiles(handles)
    global exePath;
    global exeName;

    % init
    N = length(handles.csvFiles);
    
    % process each file
    for i = 1:N
        % init data
        X = [];
        exSignal = [];
        nodeControl = [];
        exControl = [];
        
        % load node status signals csv or mat file
        fname = handles.csvFiles{i};
        if ~exist(fname,'file')
            disp(['file is not found. ignoring : ' fname]);
            continue;
        end
        [path,name,ext] = fileparts(fname);
        if strcmp(ext,'.mat')
            load(fname);
        else
            T = readtable(fname);
            X = table2array(T);
        end
        nodeNum = size(X,1);
        sigLen = size(X,2);
        
        if handles.format==2 % if save format is mat(all)
            if i==1
                savename = name;
            end
        else
            savename = name;
        end

        % signal transform raw or not
        if handles.transform == 1
            [X, sig, c, maxsi, minsi] = convert2SigmoidSignal(X, handles.transopt);
        end

        % show node status signals
        if handles.showSig > 0
            figure; plot(X.');
            title(['Node Status Signals : ' name]);
            xlabel('Time Series');
            ylabel('Signal Value');
        end

        % calc VARDNN surrogate
        if handles.vardnn > 0
            % train VARDNN network
            vardnnFile = [exePath '/results/cache-vardnn-' name '.mat'];
            if exist(vardnnFile, 'file') && handles.noCache == 0
                disp(['read cache file : ' vardnnFile]);
                load(vardnnFile);
            else
                % layer parameters
                net = initMvarDnnNetwork(X, exSignal, nodeControl, exControl, handles.lag);
                % training VARDNN network
                miniBatchSize = ceil(sigLen / 3);
                options = trainingOptions('adam', ...
                    'ExecutionEnvironment','cpu', ...
                    'MaxEpochs', handles.maxEpochs, ...
                    'MiniBatchSize', miniBatchSize, ...
                    'Shuffle', 'every-epoch', ...
                    'GradientThreshold', 5,...
                    'L2Regularization', handles.L2Regularization, ...
                    'Verbose',false);

                disp('start training');
                net = trainMvarDnnNetwork(X, exSignal, nodeControl, exControl, net, options);
                [time, loss, rsme] = getMvarDnnTrainingResult(net);
                disp(['VARDNN training result : rsme=' num2str(rsme)]);
                if handles.noCache == 0
                    save(vardnnFile, 'net');
                end
            end

            % output result matrix files
            if handles.multi > 0
                Y = surrogateMvarDnn(X, exSignal, nodeControl, exControl, net, handles.noiseType, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_vardnn_multi']);
            end
            if handles.uni > 0
                disp('vardnn_uni combination is not supported.');
            end
        end

        % calc VAR surrogate
        if handles.var > 0
            if handles.multi > 0
                net = initMvarNetwork(X, exSignal, nodeControl, exControl, handles.lag);
                Y = surrogateMVAR(X, exSignal, nodeControl, exControl, net, handles.noiseType, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_var_multi']);
            end
            if handles.uni > 0
                Y = surrogateAR(X, handles.lag, handles.noiseType, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_ar_uni']);
            end
        end

        % calc PCVAR surrogate
        if handles.pcvar > 0
            if handles.multi > 0
                net = initMpcvarNetwork(X, exSignal, nodeControl, exControl, handles.lag);
                Y = surrogateMpcvar(X, exSignal, nodeControl, exControl, net, handles.noiseType, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_pcvar_multi']);
            end
            if handles.uni > 0
                disp('pcvar_uni combination is not supported.');
            end
        end

        % calc LL surrogate
        if handles.lazy > 0
            if handles.multi > 0
                LL = initLazyLearning(X, exSignal, nodeControl, exControl, handles.lag);
                Y = surrogateLazyLearning(X, exSignal, nodeControl, exControl, LL, handles.nn, handles.noiseType, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_lazy_multi']);
            end
            if handles.uni > 0
                disp('lazy_uni combination is not supported.');
            end
        end

        % calc Random Gaussian surrogate
        if handles.rg > 0
            if handles.multi > 0
                Y = surrogateMRG(X, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_rg_multi']);
            end
            if handles.uni > 0
                Y = surrogateRG(X, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_rg_uni']);
            end
        end

        % calc Random Shuffling surrogate
        if handles.rs > 0
            if handles.multi > 0
                Y = surrogateMRS(X, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_rs_multi']);
            end
            if handles.uni > 0
                Y = surrogateRS(X, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_rs_uni']);
            end
        end

        % calc Fourier Transform surrogate
        if handles.ft > 0
            if handles.multi > 0
                Y = surrogateMFT(X, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_ft_multi']);
            end
            if handles.uni > 0
                Y = surrogateFT(X, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_ft_uni']);
            end
        end

        % calc AAFT surrogate
        if handles.aaft > 0
            if handles.multi > 0
                Y = surrogateMAAFT(X, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_aaft_multi']);
            end
            if handles.uni > 0
                Y = surrogateAAFT(X, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_aaft_uni']);
            end
        end

        % calc IAAFT surrogate
        if handles.iaaft > 0
            if handles.multi > 0
                Y = surrogateMIAAFT(X, 100, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_iaaft_multi']);
            end
            if handles.uni > 0
                Y = surrogateIAAFT(X, 100, handles.surrNum);
                saveResultFiles(handles, Y, [savename '_iaaft_uni']);
            end
        end
    end
end

%%
% output result matrix files
%
function saveResultFiles(handles, X, outname)
    if handles.format == 1
        for i=1:size(X,3)
            Y = squeeze(X(:,:,i));
            outfname = [handles.outpath '/' outname '_' num2str(i) '.mat'];
            save(outfname,'Y');
            disp(['output mat file : ' outfname]);
        end
    elseif handles.format == 2
        Y = X;
        outfname = [handles.outpath '/' outname '_all.mat'];
        save(outfname,'Y');
        disp(['output mat file : ' outfname]);
    else
        % output result matrix csv file
        for i=1:size(X,3)
            Y = squeeze(X(:,:,i));
            outputCsvFile(Y, [handles.outpath '/' outname '_' num2str(i) '.csv']);
        end
    end

    % show first sample of node status signals
    if handles.showSig > 0
        Y = squeeze(X(:,:,1));
        figure; plot(Y.');
        title(['First sample of node signals : ' outname]);
        xlabel('Time Series');
        ylabel('Signal Value');
    end
end

%%
% output csv file function
%
function outputCsvFile(mat, outfname)
    T = array2table(mat);
    writetable(T,outfname,'WriteVariableNames',false);
    disp(['output csv file : ' outfname]);
end

%%
% GSDGM command line tool

function gsdgm(varargin)

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
    handles.format = 1;
    handles.transform = 0;
    handles.transopt = NaN;
    handles.range = 'auto';
    handles.showInput = 0;
    handles.showInputRas = 0;
    handles.showSig = 0;
    handles.showRas = 0;
    
    handles.var = 0;
    handles.pcvar = 0;
    handles.vardnn = 0;

    handles.lag = 3;
    handles.noiseType = 'gaussian';
    handles.surrNum = 1;
    handles.sigLen = 0;
    handles.pcRate = 0.99;
    handles.maxEpochs = 1000;

    % load command line input
    i = 1;
    while true
        if i > size(varargin, 2)
            break;
        end
        switch varargin{i}
            case {'-v','--var'}
                handles.var = 1;
            case {'-p','--pcvar'}
                handles.pcvar = 1;
            case {'-d','--vardnn'}
                handles.vardnn = 1;
            case {'--noise'}
                handles.noiseType = varargin{i+1};
                i = i + 1;
            case {'--surrnum'}
                handles.surrNum = str2num(varargin{i+1});
                i = i + 1;
            case {'--siglen'}
                handles.sigLen = str2num(varargin{i+1});
                i = i + 1;
            case {'--lag'}
                handles.lag = str2num(varargin{i+1});
                i = i + 1;
            case {'--pcrate'}
                handles.pcRate = str2num(varargin{i+1});
                i = i + 1;
            case {'--epoch'}
                handles.maxEpochs = str2num(varargin{i+1});
                i = i + 1;
            case {'--outpath'}
                handles.outpath = varargin{i+1};
                i = i + 1;
            case {'--format'}
                handles.format = str2num(varargin{i+1});
                i = i + 1;
            case {'--range'}
                handles.range = varargin{i+1};
                i = i + 1;
            case {'--transform'}
                handles.transform = str2num(varargin{i+1});
                i = i + 1;
            case {'--transopt'}
                handles.transopt = str2num(varargin{i+1});
                i = i + 1;
            case {'--showinsig'}
                handles.showInput = 1;
            case {'--showinras'}
                handles.showInputRas = 1;
            case {'--showsig'}
                handles.showSig = 1;
            case {'--showras'}
                handles.showRas = 1;
            case {'-h','--help'}
                showUsage();
                return;
            case {'--version'}
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
        disp('no input files. please specify time-series files.');
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
    disp(['model training : ' exeName ' [options] file1.mat file2.mat ...']);
    disp(['surrogate data : ' exeName ' [options] file_gsm_<type>.mat']);
    disp('  -v, --var           output Vector Auto-Regression (VAR) group surrogate model (<filename>_gsm_var.mat)');
    disp('  -p, --pcvar         output Principal Component VAR (PCVAR) group surrogate model (<filename>_gsm_pcvar.mat)');
    disp('  -d, --vardnn        output VAR Deep Neural Network (VARDNN) group surrogate model (<filename>_gsm_vardnn.mat)');
    disp('  --lag num           time lag <num> for VAR, PCVAR, VARDNN surrogate model (default:3)');
    disp('  --noise type        noise type for VAR, PCVAR, VARDNN surrogate model (default:"gaussian" or "residuals")');
    disp('  --outpath path      output files <path> (default:"results")');
    disp('  --transform type    input training signal transform <type> 0:raw, 1:sigmoid (default:0)');
    disp('  --transopt num      signal transform option <num> (for type 1:centroid value)');
    disp('  --format type       output surrogate data file format <type> 0:csv, 1:mat (default:1)');
    disp('  --surrnum num       output surrogate sample number <num> (default:1)');
    disp('  --siglen num        output time-series length <num> (default:same as input time-series)');
    disp('  --range type        output surrogate value range (default:"auto", sigma:<num>, full:<num>, <min>:<max> or "none")');
    disp('  --pcrate num        principal component variance rate <num> for PCVAR surrogate (default:0.99)');
    disp('  --epoch num         VARDNN surrogate training epoch number <num> (default:1000)');
    disp('  --showinsig         show input time-series data of <filename>.csv');
    disp('  --showinras         show raster plot of input time-series data of <filename>.csv');
    disp('  --showsig           show output surrogate time-series data');
    disp('  --showras           show raster plot of output surrogate time-series data');
    disp('  --version           show version number');
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

    % load each file
    CX = {}; names = {}; net = []; gRange = []; savename = '';
    for i = 1:N
        argv = handles.csvFiles{i};
        % check url or file
        if contains(argv, 'http://') || contains(argv, 'https://')
            % make download cache directory
            if ~exist('data/cache','dir')
                mkdir('data/cache');
            end
            % make cache file string
            url = argv;
            argv = ['data/cache/' url2cacheString(argv)];
            if ~exist(argv,'file')
                disp(['downloading ' url ' ...']);
                websave(argv, url);
                disp(['save cache file : ' argv]);
            end
        end
        
        % load multivariate time-series csv or mat file
        flist = dir(argv);
        if isempty(flist)
            disp(['file is not found. ignoring : ' argv]);
            continue;
        end
        for k=1:length(flist)
            % init data
            X = [];

            fname = [flist(k).folder '/' flist(k).name];
            [path,name,ext] = fileparts(fname);
            if strcmp(ext,'.mat')
                f = load(fname);
                if isfield(f,'CX')
                    % training mode
                    if isfield(f,'multiple') && isa(f.CX{1},'uint16')
                        % uint16 for demo
                        tn = cell(1,length(f.CX));
                        for j=1:length(f.CX)
                            tn{j} = single(f.CX{j}) / f.multiple;
                        end
                        CX = [CX, tn];
                    else
                        CX = [CX, f.CX]; % single
                    end
                    if isfield(f,'names')
                        tn = cell(1,length(f.CX));
                        for j=1:length(f.CX)
                            tn{j} = strrep(f.names{j},'_','-');
                        end
                        names = [names, tn];
                    else
                        tn = {};
                        for j=1:length(f.CX)
                            tn{j} = [strrep(name,'_','-') '-' num2str(j)];
                        end
                        names = [names, tn];
                    end
                elseif isfield(f,'X')
                    % training mode
                    if isfield(f,'name'), name = f.name; end
                    names = [names, strrep(name,'_','-')];
                    CX = [CX, f.X];
                elseif isfield(f,'net')
                    % surrogate data mode
                    if isfield(f,'name'), name = f.name; end
                    if isfield(f,'gRange'), gRange = f.gRange; end
                    net = f.net;
                else
                    disp(['file does not contain "X" matrix or "CX" cell. ignoring : ' fname]);
                end
            else
                % training mode
                T = readtable(fname);
                X = table2array(T);
                names = [names, strrep(name,'_','-')];
                CX = [CX, X];
            end

            if isempty(savename)
                savename = name;
            end
        end
    end
    
    % check each multivariate time-series
    for i = 1:length(CX)
        X = CX{i};
        
        % signal transform raw or not
        if handles.transform == 1
            [X, sig, c, maxsi, minsi] = convert2SigmoidSignal(X, handles.transopt);
            CX{i} = X;
        end
        
        % show input signals
        if handles.showInput > 0
            figure; plot(X.');
            title(['Input Signals : ' names{i}]);
            xlabel('Time Series');
            ylabel('Signal Value');
        end
        % show input signals
        if handles.showInputRas > 0
            figure; imagesc(X);
            title(['Raster plot of Input signals : ' names{i}]);
            xlabel('Time Series');
            ylabel('Node number');
            colorbar;
        end
    end
    
    % training mode
    if ~isempty(CX)
        % get group range
        gRange = getGroupRange(CX);

        if handles.var > 0
            net = initMvarNetworkWithCell(CX, [], [], [], handles.lag);
            saveModelFile(handles, net, gRange, [savename '_gsm_var']);
        end
        if handles.pcvar > 0
            net = initMpcvarNetworkWithCell(CX, [], [], [], handles.lag, handles.pcRate);    
            saveModelFile(handles, net, gRange, [savename '_gsm_pcvar']);
        end
        if handles.vardnn > 0
            % set training options
            sigLen = size(CX{1},2);
            miniBatchSize = ceil(sigLen / 3);

            options = trainingOptions('adam', ...
                'ExecutionEnvironment','cpu', ...
                'MaxEpochs',handles.maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'Shuffle','every-epoch', ...
                'GradientThreshold',5,...
                'Verbose',false);

            net = initMvarDnnNetworkWithCell(CX, [], [], [], handles.lag, 60);
            net = trainMvarDnnNetworkWithCell(CX, [], [], [], net, options);
            saveModelFile(handles, net, gRange, [savename '_gsm_vardnn']);
        end
    end
    
    % surrogate data mode
    if ~isempty(net)
        if handles.sigLen > 0, sigLen = handles.sigLen; else sigLen = net.sigLen; end

        % dummy signal for nodeNum, sigLen and surrogate initial value (also affect surrogate value range)
        X = (mvnrnd(net.cxM, net.cxCov, sigLen))';
        
        % set output value range
        range = NaN; % unknown. calc range based on X.
        if strcmp(handles.range,'auto')
            % 3.6 sigma of the whole group
            if ~isempty(gRange)
                range = [gRange.m - gRange.s * 3.6, gRange.m + gRange.s * 3.6];
            end
        elseif strcmp(handles.range,'none')
            range = []; % empty. no range limit
        elseif contains(handles.range,':')
            str = split(handles.range,':');
            if strcmp(str{1},'sigma') % <num> sigma of the whole group
                if ~isempty(gRange)
                    n = str2num(str{2});
                    range = [gRange.m - gRange.s * n, gRange.m + gRange.s * n];
                end
            elseif strcmp(str{1},'full') % <num> * full min & max range of the whole group
                if ~isempty(gRange)
                    n = (str2num(str{2}) - 1) / 2;
                    r = gRange.max - gRange.min;
                    range = [-gRange.min - r*n, gRange.max + r*n];
                end
            else
                % force [<num>, <num>] range
                range = [str2num(str{1}),str2num(str{2})];
            end
        end
        
        % generate surrogate data
        if isfield(net,'nodeNetwork')
            nettype = 'vardnn';
            Y = surrogateMvarDnn(X, [], [], [], net, handles.noiseType, handles.surrNum, range);
        elseif isfield(net,'mu')
            nettype = 'pcvar';
            Y = surrogateMpcvar(X, [], [], [], net, handles.noiseType, handles.surrNum, range);
        else
            nettype = 'var';
            Y = surrogateMVAR(X, [], [], [], net, handles.noiseType, handles.surrNum, range);
        end

        CX = cell(1,size(Y,3)); names = cell(1,size(Y,3));
        for i=1:length(CX)
            CX{i} = squeeze(Y(:,:,i));
            names{i} = [savename '-gsd-' nettype '-' num2str(i)];

            % show output signals
            if handles.showSig > 0
                figure; plot(CX{i}.');
                title(['Group Surrogate Data : ' strrep(names{i},'_','-')]);
                xlabel('Time Series');
                ylabel('Signal Value');
            end

            % show output signals
            if handles.showRas > 0
                figure; imagesc(CX{i});
                title(['Raster plot of Group Surrogate Data : ' strrep(names{i},'_','-')]);
                xlabel('Time Series');
                ylabel('Node number');
                colorbar;
            end
        end

        % output result matrix files
        sn = strrep(savename,'_gsm_vardnn','');
        sn = strrep(sn,'_gsm_var','');
        sn = strrep(sn,'_gsm_pcvar','');
        saveResultFiles(handles, CX, names, [sn '_gsd_' nettype]);
    end
end

%%
% output result files
%
function saveModelFile(handles, net, gRange, outname)
    outfname = [handles.outpath '/' outname '.mat'];
    save(outfname, 'net', 'gRange', '-v7.3');
    disp(['output group surrogate model file : ' outfname]);
end

function saveResultFiles(handles, CX, names, outname)
    if handles.format == 1
        outfname = [handles.outpath '/' outname '.mat'];
        save(outfname, 'CX', 'names', '-v7.3');
        disp(['output mat file : ' outfname]);
    else
        % output result matrix csv file
        for i=1:length(CX)
            outputCsvFile(CX{i}, [handles.outpath '/' outname '_' num2str(i) '.csv']);
        end
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

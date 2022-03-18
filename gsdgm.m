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
    handles.showInput = 0;
    handles.showSig = 0;
    
    handles.var = 0;
    handles.pcvar = 0;
    handles.vardnn = 0;

    handles.lag = 3;
    handles.noiseType = 'gaussian';
    handles.pcRate = 0.99;
    handles.surrNum = 1;
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
            case {'--transform'}
                handles.transform = str2num(varargin{i+1});
                i = i + 1;
            case {'--transopt'}
                handles.transopt = str2num(varargin{i+1});
                i = i + 1;
            case {'--showinsig'}
                handles.showInput = 1;
            case {'--showsig'}
                handles.showSig = 1;
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
    disp('  --pcrate num        principal component variance rate <num> for PCVAR surrogate (default:0.99)');
    disp('  --epoch num         VARDNN surrogate training epoch number <num> (default:1000)');
    disp('  --showinsig         show each time-series data of <filename>.csv');
    disp('  --showsig           show output surrogate time-series data');
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
    CX = {}; names = {}; net = [];
    for i = 1:N
        % init data
        X = [];

        % load node status signals csv or mat file
        fname = handles.csvFiles{i};
        if ~exist(fname,'file')
            disp(['file is not found. ignoring : ' fname]);
            continue;
        end
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
                names = [names, strrep(name,'_','-')];
                CX = [CX, f.X];
            elseif isfield(f,'net')
                % surrogate data mode
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
        
        if i==1
            savename = name;
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
            title(['Input signals : ' names{i}]);
            xlabel('Time Series');
            ylabel('Signal Value');
        end
    end
    
    % training mode
    if ~isempty(CX)
        if handles.var > 0
            net = initMvarNetworkWithCell(CX, [], [], [], handles.lag);
            saveModelFile(handles, net, [savename '_gsm_var']);
        end
        if handles.pcvar > 0
            net = initMpcvarNetworkWithCell(CX, [], [], [], handles.lag, handles.pcRate);    
            saveModelFile(handles, net, [savename '_gsm_pcvar']);
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
            saveModelFile(handles, net, [savename '_gsm_vardnn']);
        end
    end
    
    % surrogate data mode
    if ~isempty(net)
        % dummy signal for nodeNum, sigLen and surrogate initial value (also affect surrogate value range)
        X = (mvnrnd(net.cxM, net.cxCov, net.sigLen))';

        % generate surrogate data
        if isfield(net,'nodeNetwork')
            nettype = 'vardnn';
            Y = surrogateMvarDnn(X, [], [], [], net, handles.noiseType, handles.surrNum);
        elseif isfield(net,'mu')
            nettype = 'pcvar';
            Y = surrogateMpcvar(X, [], [], [], net, handles.noiseType, handles.surrNum);
        else
            nettype = 'var';
            Y = surrogateMVAR(X, [], [], [], net, handles.noiseType, handles.surrNum);
        end

        CX = cell(1,size(Y,3)); names = cell(1,size(Y,3));
        for i=1:length(CX)
            CX{i} = squeeze(Y(:,:,i));
            names{i} = [savename '-gsd-' nettype '-' num2str(i)];
        
            % show input signals
            if handles.showSig > 0
                figure; plot(CX{i}.');
                title(['Group Surrogate Data : ' names{i}]);
                xlabel('Time Series');
                ylabel('Signal Value');
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
function saveModelFile(handles, net, outname)
    outfname = [handles.outpath '/' outname '.mat'];
    save(outfname, 'net', '-v7.3');
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

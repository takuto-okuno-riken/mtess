%%
% MTESS command line tool

function mtess(varargin)

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
    handles.range = NaN;
    handles.nDft = 100;
    handles.pcc = 0;
    handles.cclag = 8;
    handles.pcclag = 8;
    handles.outpath = 'results';
    handles.format = 0;
    handles.transform = 0;
    handles.transopt = NaN;
    handles.showInput = 0;
    handles.showMat = 0;
    handles.showSig = 0;
    handles.showProp = 0;
    handles.showNode = 0;
    handles.cache = '';

    % load command line input
    i = 1;
    while true
        if i > size(varargin, 2)
            break;
        end
        switch varargin{i}
            case {'--range'}
                str = strsplit(varargin{i+1},':');
                handles.range = [str2num(str{1}) str2num(str{2})];
                i = i + 1;
            case {'--ndft'}
                handles.nDft = str2num(varargin{i+1});
                i = i + 1;
            case {'--pcc'}
                handles.pcc = str2num(varargin{i+1});
                i = i + 1;
            case {'--cclag'}
                handles.cclag = str2num(varargin{i+1});
                i = i + 1;
            case {'--pcclag'}
                handles.pcclag = str2num(varargin{i+1});
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
            case {'--showmat'}
                handles.showMat = 1;
            case {'--showprop'}
                handles.showProp = 1;
            case {'--shownode'}
                handles.showNode = 1;
            case {'--cache'}
                handles.cache = varargin{i+1};
                i = i + 1;
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
    disp(['usage: ' exeName ' [options] file1.csv file2.csv ...']);
    disp('  --range n1:n2       value range [n1, n2] for normalized mean and std dev (default:min and max of input data)');
    disp('  --ndft num          DFT sampling <number> (even number) (default: 100)');
    disp('  --pcc type          Partial Cross-Correlation algorithm 0:auto, 1:PCC, 2:SV-PCC (dafault:0)');
    disp('  --cclag num         time lag <num> for Cross Correlation (default:8)');
    disp('  --pcclag num        time lag <num> for Partial Cross Correlation (default:8)');
    disp('  --outpath           output files path (default:"results")');
    disp('  --format type       save file format <type> 0:csv, 1:mat (default:0)');
    disp('  --transform type    input signal transform <type> 0:raw, 1:sigmoid (default:0)');
    disp('  --transopt num      signal transform option <num> (for type 1:centroid value)');
    disp('  --showinsig         show input signals of <filename>.csv');
    disp('  --showmat           show result MTESS matrix');
    disp('  --showsig           show 1 vs. others node signals');
    disp('  --showprop          show result polar chart of 1 vs. others MTESS statistical properties');
    disp('  --shownode          show result line plot of 1 vs. others node MTESS');
    disp('  --cache filename    use cache <filename> for MTESS calculation');
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

    % load each file
    CX = {}; names = {}; issetname = false;
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
                CX = f.CX;
                if isfield(f,'names')
                    names = f.names;
                end
            elseif isfield(f,'X')
                names{i} = strrep(name,'_','-');
                CX{i} = f.X;
            else
                disp(['file does not contain "X" matrix or "CX" cell. ignoring : ' fname]);
            end
        else
            T = readtable(fname);
            X = table2array(T);
            names{i} = strrep(name,'_','-');
            CX{i} = X;
        end
        
        if i==1
            savename = name;
        end
    end

    % check name
    if isempty(names)
        issetname = true;
    end
    
    % check each multivariate time-series
    nodeNum = size(CX{1},1);
    xmax = NaN;
    xmin = NaN;
    for i = 1:length(CX)
        X = CX{i};
        
        % signal transform raw or not
        if handles.transform == 1
            [X, sig, c, maxsi, minsi] = convert2SigmoidSignal(X, handles.transopt);
            CX{i} = X;
        end

        if isnan(xmax) || xmax < max(X,[],'all')
            xmax = max(X,[],'all');
        end
        if isnan(xmin) || xmin > min(X,[],'all')
            xmin = min(X,[],'all');
        end
        
        % check name
        if issetname
            names{i} = ['data' num2str(i)];
        end

        % show input signals
        if handles.showInput > 0
            figure; plot(X.');
            title(['File Signals : ' names{i}]);
            xlabel('Time Series');
            ylabel('Signal Value');
        end
    end
    if isnan(handles.range)
        handles.range = [xmin, xmax];
    end
    
    % calc MTESS
    if handles.pcc == 1
        pccFunc = @calcPartialCrossCorrelation;
    elseif handles.pcc == 2
        pccFunc = @calcSvPartialCrossCorrelation;
    else
        if nodeNum < 48
            pccFunc = @calcPartialCrossCorrelation;
        else
            pccFunc = @calcSvPartialCrossCorrelation;
        end
    end
    [MTS, MTSp, nMTS, nMTSp, Means, Stds, Amps, FCs, PCs, CCs, PCCs] = calcMtess(CX, handles.range, handles.nDft, pccFunc, handles.cclag, handles.pcclag, handles.cache);

    % output result matrix files
    saveResultFiles(handles, MTS, MTSp, nMTS, nMTSp, savename);

    % show all matrix
    if handles.showMat > 0
        plotMtessAllMatrix(MTS, MTSp, savename);
    end
    
    % show 1 vs. others signals
    if handles.showSig > 0
        for i = 2:length(CX)
            figure; plotTwoSignals(CX{1},CX{i},0,handles.range);
            sgt = sgtitle(['Node signals : ' names{1} ' vs. ' names{i}]);
            sgt.FontSize = 10;
            legend(names([1,i]));
        end
    end
    
    % show 1 vs. others MTESS statistical properties
    if handles.showProp > 0
        P=squeeze(MTSp(1,2:length(CX),:));
        figure; plotMtessSpiderPlot(P);
        legend(names(2:length(CX)));
        title('MTESS polar chart : 1 vs. others');
    end
    
    % show 1 vs. others node MTESS
    if handles.showNode > 0
        A = squeeze(nMTS(1,2:length(CX),:));
        figure; plot(A.'); title('node MTESS : 1 vs. others');
        legend(names(2:length(CX))); ylim([0,5]);
        xlabel('Node number');
        ylabel('MTESS');
    end
end

%%
% output result matrix files
%
function saveResultFiles(handles, MTS, MTSp, nMTS, nMTSp, outname)
    if handles.format == 1
        save([handles.outpath '/' outname '_mtess.mat'], 'MTS', 'MTSp', 'nMTS', 'nMTSp');
    else
        % output result MTESS matrix csv file
        outputCsvFile(MTS, [handles.outpath '/' outname '_mtess.csv']);

        % output result MTESS statistical property matrix csv file
        props = {'M','SD','Amp','FC','PC','CC','PCC'};
        for i=1:length(props)
            outputCsvFile(MTSp(:,:,i), [handles.outpath '/' outname '_mtess_' props{i} '.csv']);
        end

        % output result node MTESS matrix csv file
        props = {'M','SD','Amp','FC','PC','CC','PCC'};
        for i=1:size(nMTS,3)
            outputCsvFile(nMTS(:,:,i), [handles.outpath '/' outname '_mtess_node' num2str(i) '.csv']);
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

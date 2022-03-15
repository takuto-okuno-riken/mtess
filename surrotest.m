%%
% surrogate test command line tool

function surrotest(varargin)

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
    handles.format = 0;
    handles.showSig = 0;
    handles.showRank = 0;

    handles.gaussian = 0;
    handles.linear = 0;
    handles.iid = 0;
    handles.side = 2;

    % load command line input
    i = 1;
    while true
        if i > size(varargin, 2)
            break;
        end
        switch varargin{i}
            case {'-g','--gaussian'}
                handles.gaussian = 1;
            case {'-l','--linear'}
                handles.linear = 1;
            case {'-i','--iid'}
                handles.iid = 1;
            case {'--outpath'}
                handles.outpath = varargin{i+1};
                i = i + 1;
            case {'--format'}
                handles.format = str2num(varargin{i+1});
                i = i + 1;
            case {'--side'}
                handles.side = str2num(varargin{i+1});
                i = i + 1;
            case {'--showsig'}
                handles.showSig = 1;
            case {'--showrank'}
                handles.showRank = 1;
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
        disp('original and surrogate files are required. please specify node status signal files.');
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
    
    disp(['usage: ' exeName ' [options] <original>.csv surrogate.mat ...']);
    disp('  -g, --gaussian      output Gaussian distribution test (<original>_gauss_test.csv)');
    disp('  -l, --linear        output Linearity test  (<original>_linear_test.csv)');
    disp('  -i, --iid           output I.I.D test (<original>_iid_test.csv)');
    disp('  --side num          bottm-side(1), both-side(2), top-side(3) (default:2)');
    disp('  --outpath           output files path (default:"results")');
    disp('  --format type       save file format <type> 0:csv, 1:mat (default:0)');
    disp('  --showsig           show node status signals of <original>.csv');
    disp('  --showrank          show rank result of <original>.csv');
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
    X = []; % original
    Y = []; % surrogate
    
    % process each file
    for i = 1:N
        % load node status signals csv or mat file
        fname = handles.csvFiles{i};
        if ~exist(fname,'file')
            disp(['file is not found. ignoring : ' fname]);
            continue;
        end
        [path,name,ext] = fileparts(fname);
        if strcmp(ext,'.mat')
            f = load(fname);
            if i==1
                X = f.X;
                if isfield(f,'Y')
                    Y = f.Y;
                end
            else
                Y = cat(3,Y,f.Y);
            end
        else
            if i==1
                T = readtable(fname);
                X = table2array(T);
            else
                T = readtable(fname);
                Y1 = table2array(T);
                Y = cat(3,Y,Y1);
            end
        end

        if i==1
            savename = name;
        end
    end

    requireNum = 39;
    if handles.side == 1 || handles.side == 3
        requireNum = 19;
    end
    if size(Y,3) < requireNum
        disp(['The surrogate test requires at least ' num2str(requireNum) ' surrogate data.']);
        return;
    end
    nodeNum = size(X,1);

    % show node status signals
    if handles.showSig > 0
        figure; plot(X.');
        title(['Node Status Signals : ' name]);
        xlabel('Time Series');
        ylabel('Signal Value');
    end

    % surrogate linearity test
    if handles.linear > 0
        statisticParams = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
        if handles.showRank > 0
            for i=1:nodeNum
                figure; [H(i), P(i), ~, Rank(i)] = plotSurrogateTest(X(i,:), Y(i,:,:), @calcSurrStatLinear, statisticParams, handles.side);
                title(['Linear surrogate test Node ' num2str(i) ' : rank=' num2str(Rank(i)) ', p-value=' num2str(P(i))]);
            end
            P=P.'; Rank=Rank.';
        else
            [H, P, ~, Rank] = calcSurrogateTest(X, Y, @calcSurrStatLinear, statisticParams, handles.side);
        end
        disp(['significantly not linear (' num2str(sum(H)) ' / ' num2str(nodeNum) ')']);
        saveResultFiles(handles, P, Rank, [savename '_linear_test']);
    end

    % surrogate gaussian distribution test
    if handles.gaussian > 0
        if handles.showRank > 0
            for i=1:nodeNum
                figure; [H(i), P(i), ~, Rank(i)] = plotSurrogateTest(X(i,:), Y(i,:,:), @calcSurrStatGaussian, [], handles.side);
                title(['Gaussian distribution surrogate test Node ' num2str(i) ' : rank=' num2str(Rank(i)) ', p-value=' num2str(P(i))]);
            end
            P=P.'; Rank=Rank.';
        else
            [H, P, ~, Rank] = calcSurrogateTest(X, Y, @calcSurrStatGaussian, [], handles.side);
        end
        disp(['significantly not gaussian (' num2str(sum(H)) ' / ' num2str(nodeNum) ')']);
        saveResultFiles(handles, P, Rank, [savename '_gaussian_test']);
    end
    
    % surrogate I.I.D test
    if handles.iid > 0
        if handles.showRank > 0
            for i=1:nodeNum
                figure; [H(i), P(i), ~, Rank(i)] = plotSurrogateTest(X(i,:), Y(i,:,:), @calcSurrStatIID, [], handles.side);
                title(['I.I.D surrogate test Node ' num2str(i) ' : rank=' num2str(Rank(i)) ', p-value=' num2str(P(i))]);
            end
            P=P.'; Rank=Rank.';
        else
            [H, P, ~, Rank] = calcSurrogateTest(X, Y, @calcSurrStatIID, [], handles.side);
        end
        disp(['significantly not I.I.D (' num2str(sum(H)) ' / ' num2str(nodeNum) ')']);
        saveResultFiles(handles, P, Rank, [savename '_iid_test']);
    end
end

%%
% output result matrix files
%
function saveResultFiles(handles, P, Rank, outname)
    if handles.format == 1
        outfname = [handles.outpath '/' outname '.mat'];
        save(outfname,'P', 'Rank');
        disp(['output mat file : ' outfname]);
    else
        % output result matrix csv file
        outputCsvFile(P, [handles.outpath '/' outname '_pval.csv']);
        outputCsvFile(Rank, [handles.outpath '/' outname '_rank.csv']);
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

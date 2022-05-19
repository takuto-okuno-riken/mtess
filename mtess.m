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
    handles.range = 'auto';
    handles.nDft = 100;
    handles.pcc = 0;
    handles.cclag = NaN;
    handles.pcclag = NaN;
    handles.outpath = 'results';
    handles.format = 1;
    handles.transform = 0;
    handles.transopt = NaN;
    handles.showInput = 0;
    handles.showInputRas = 0;
    handles.showMat = 0;
    handles.showSig = 0;
    handles.showProp = 0;
    handles.showNode = 0;
    handles.showDend = '';
    handles.showForce = 0;
    handles.noCache = 0;

    % load command line input
    i = 1;
    while true
        if i > size(varargin, 2)
            break;
        end
        switch varargin{i}
            case {'--range'}
                handles.range = varargin{i+1};
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
            case {'--showinras'}
                handles.showInputRas = 1;
            case {'--showsig'}
                handles.showSig = 1;
            case {'--showmat'}
                handles.showMat = 1;
            case {'--showprop'}
                handles.showProp = 1;
            case {'--shownode'}
                handles.showNode = 1;
            case {'--showforce'}
                handles.showForce = 1;
            case {'--showdend'}
                handles.showDend = varargin{i+1};
                i = i + 1;
            case {'--nocache'}
                handles.noCache = 1;
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
    disp(['usage: ' exeName ' [options] file1.mat file2.mat ...']);
    disp('  --range type        input group value range (default:"auto", sigma:<num>, full:<num> or <min>:<max>)');
    disp('  --ndft num          DFT sampling <number> (even number) (default: 100)');
    disp('  --pcc type          Partial Cross-Correlation algorithm 0:auto, 1:PCC, 2:SV-PCC, 3:PC-PCC (dafault:0)');
    disp('  --cclag num         time lag <num> for Cross Correlation (default:8)');
    disp('  --pcclag num        time lag <num> for Partial Cross Correlation (default:8)');
    disp('  --outpath path      output files <path> (default:"results")');
    disp('  --format type       save file format <type> 0:csv, 1:mat (default:1)');
    disp('  --transform type    input signal transform <type> 0:raw, 1:sigmoid (default:0)');
    disp('  --transopt num      signal transform option <num> (for type 1:centroid value)');
    disp('  --showinsig         show input time-series data of <filename>.csv');
    disp('  --showinras         show raster plot of input time-series data of <filename>.csv');
    disp('  --showmat           show result MTESS matrix');
    disp('  --showsig           show 1 vs. others node signals');
    disp('  --showprop          show result polar chart of 1 vs. others MTESS statistical properties');
    disp('  --shownode          show result line plot of 1 vs. others node MTESS');
    disp('  --showdend algo     show dendrogram of <algo> hierarchical clustering based on MTESS matrix. see MATLAB linkage method option.');
    disp('  --showforce         show force weight effect graph based on MTESS matrix');
    disp('  --nocache           do not use cache file for MTESS calculation');
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
    CX = {}; names = {}; savename = '';
    for i = 1:N
        % load multivariate time-series csv or mat file
        argv = handles.csvFiles{i};
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
                        tn = cell(1,length(f.CX));
                        for j=1:length(f.CX)
                            tn{j} = [strrep(name,'_','-') '-' num2str(j)];
                        end
                        names = [names, tn];
                    end
                elseif isfield(f,'X')
                    names = [names, strrep(name,'_','-')];
                    CX = [CX, f.X];
                else
                    disp(['file does not contain "X" matrix or "CX" cell. ignoring : ' fname]);
                end
            else
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
    nodeNum = size(CX{1},1);
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
            title(['Raster plot of Signals : ' names{i}]);
            xlabel('Time Series');
            ylabel('Node number');
            colorbar;
        end
    end

    % get group range
    gRange = getGroupRange(CX);

    % set group value range
    range = NaN; % unknown. calc range based on X.
    if strcmp(handles.range,'auto')
        % 3 sigma of the whole group
        range = [gRange.m - gRange.s * 3, gRange.m + gRange.s * 3];
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
                range = [gRange.min - r*n, gRange.max + r*n];
            end
        else
            % force [<num>, <num>] range
            range = [str2num(str{1}),str2num(str{2})];
        end
    else
        disp('bad range option. stop operation.');
        showUsage();
        return;
    end
    
    % calc MTESS
    pcName = 'PC';
    pccFunc = @calcPartialCrossCorrelation;
    cclag = 8;
    pcclag = 8;
    if handles.pcc == 1
        % same as default
    elseif handles.pcc == 2
        pccFunc = @calcSvPartialCrossCorrelation;
        pcName = 'SVgPC';
        pcclag = 2;
    elseif handles.pcc == 3
        pccFunc = @calcPcPartialCrossCorrelation;
        pcName = 'PCPC';
    else
        % auto
        if nodeNum > 400
            % too big. use normal PCC
            pcclag = 2;
        elseif nodeNum >= 48
            pccFunc = @calcSvPartialCrossCorrelation;
            pcName = 'SVgPC';
            pcclag = 2;
        end
    end
    if ~isnan(handles.cclag), cclag = handles.cclag; end
    if ~isnan(handles.pcclag), pcclag = handles.pcclag; end

    if handles.noCache > 0
        cache = {};
    else
        cache = names;
    end
    [MTS, MTSp, nMTS, nMTSp, Means, Stds, Amps, FCs, PCs, CCs, PCCs] = calcMtess(CX, range, handles.nDft, pccFunc, cclag, pcclag, cache);

    % output result matrix files
    saveResultFiles(handles, MTS, MTSp, nMTS, nMTSp, savename);

    % show all matrix
    if handles.showMat > 0
        figure; h=bar3(MTS); title('MTESS matrix 3D bar graph');
        for i=1:length(h) h(i).FaceAlpha=0.6; h(i).EdgeAlpha=0.6; end
        xlabel('Cell number');
        ylabel('Cell number');
        zlabel('MTESS'); zticks([0 1 2 3 4 5]);
        plotMtessAllMatrix(MTS, MTSp, savename);
    end
    
    % show 1 vs. others signals
    if handles.showSig > 0
        for i = 2:length(CX)
            figure; plotTwoSignals(single(CX{1}),single(CX{i}),0,handles.range);
            sgt = sgtitle(['Node signals : ' names{1} ' vs. ' names{i}]);
            sgt.FontSize = 10;
            legend(names([1,i]));
        end
    end
    
    % show 1 vs. others MTESS statistical properties
    if handles.showProp > 0
        P=squeeze(MTSp(1,2:length(CX),:));
        figure; plotMtessSpiderPlot(P, cclag, pcclag, pcName);
        legend(names(2:length(CX)));
        title('MTESS polar chart : 1 vs. others');
    end
    
    % show 1 vs. others node MTESS
    if handles.showNode > 0
        A = squeeze(nMTS(1,2:length(CX),:));
        figure; plot(A.',':o','LineWidth',1, 'MarkerFaceColor','auto', 'MarkerSize',4);
        title('node MTESS : 1 vs. others');
        yticks([0 1 2 3 4 5]);
        legend(names(2:length(CX))); ylim([0,5]);
        xlabel('Node number');
        ylabel('MTESS');
    end
    
    % show dendrogram
    if ~isempty(handles.showDend)
        X = 5 - MTS;
        X(isnan(X))=0; X = X + X.';
        y = squareform(X);
        Z = linkage(y,'ward');
        figure; dendrogram(Z);
        title('Hierarchical clustering based on MTESS');
        ylabel('MTESS distance');
        xlabel('Cell number');
    end
    
    % show force weight effect graph
    if handles.showForce > 0
        X = 5 - MTS;
        dn = cell(1,length(CX));
        for i=1:length(dn), dn{i}=[num2str(i)]; end
        G = graph(X, dn, 'omitselfloops','upper');
        figure; gp=plot(G,'Layout','force','WeightEffect','direct');%'force','UseGravity',true);
        gp.EdgeColor = [0.7, 0.7, 0.7];
        gp.LineStyle = ':';
        title('Force weight effect graph based on MTESS');
    end
end

%%
% output result matrix files
%
function saveResultFiles(handles, MTS, MTSp, nMTS, nMTSp, outname)
    if handles.format == 1
        outfname = [handles.outpath '/' outname '_mtess.mat'];
        save(outfname, 'MTS', 'MTSp', 'nMTS', 'nMTSp', '-v7.3');
        disp(['output mat file : ' outfname]);
    else
        % output result MTESS matrix csv file
        outputCsvFile(MTS, [handles.outpath '/' outname '_mtess.csv']);

        % output result MTESS statistical property matrix csv file
        props = {'M','SD','AC','CM','PCM','CCM','PCCM'};
        for i=1:length(props)
            outputCsvFile(MTSp(:,:,i), [handles.outpath '/' outname '_mtess_' props{i} '.csv']);
        end

        % output result node MTESS matrix csv file
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

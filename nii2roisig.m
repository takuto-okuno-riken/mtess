%%
% nifti to ROI signal conversion command line tool

function nii2roisig(varargin)

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
    handles.niiFiles = {};
    handles.atlasFile = [];
    handles.outpath = 'results';
    handles.format = 2;
    handles.transform = 0;
    handles.transopt = NaN;

    handles.showSig = 0;
    handles.showRas = 0;
    handles.noCache = 0;

    % load command line input
    i = 1;
    while true
        if i > size(varargin, 2)
            break;
        end
        switch varargin{i}
            case {'-a','--atlas'}
                handles.atlasFile = varargin{i+1};
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
            case {'--showsig'}
                handles.showSig = 1;
            case {'--showras'}
                handles.showRas = 1;
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
                    handles.niiFiles = [handles.niiFiles varargin{i}];
                end
        end
        i = i + 1;
    end
    
    % check command input
    if handles.commandError
        showUsage();
        return;
    elseif isempty(handles.niiFiles) || isempty(handles.atlasFile)
        disp('ROI atlas and nifti files are required. please specify nifti files.');
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
    
    disp(['usage: ' exeName ' [options] -a atlas.nii file1.nii ...']);
    disp('  -a, --atlas file    ROI atlas nifti <file>');
    disp('  --outpath path      output files <path> (default:"results")');
    disp('  --format type       save file format <type> 0:csv, 1:mat(each), 2:mat(all) (default:2)');
    disp('  --transform type    output signal transform <type> 0:raw, 1:sigmoid (default:0)');
    disp('  --transopt num      signal transform option <num> (for type 1:centroid value)');
    disp('  --showsig           show output time-series data of <original>.csv');
    disp('  --showras           show raster plot of output time-series data of <original>.csv');
    disp('  --nocache           do not use cache file for conversion');
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
    N = length(handles.niiFiles);
    
    % load ROI atlas file
    if ~exist(handles.atlasFile,'file')
        disp(['atlas file is not found : ' handles.atlasFile]);
        return;
    end
    [~,name,~] = fileparts(handles.atlasFile);
    atlasName = strrep(name,'.nii','');

    atlasInfo = niftiinfo(handles.atlasFile);
    atlasV = niftiread(atlasInfo);
    atlasV = adjustVolumeDir(atlasV, atlasInfo);
    R = unique(atlasV);
    R(R==0) = []; % remove 0
    nodeNum = length(R);
    
    % check and make cache directory
    if ~exist('results/cache','dir')
        mkdir('results/cache');
    end

    % atlas resampling
    disp(['checking atlas space size ...']);
    name = '';
    for i = 1:N
        argv = handles.niiFiles{i};
        flist = dir(argv);
        for k=1:length(flist)
            % load multivariate time-series csv or mat file
            fname = [flist(k).folder '/' flist(k).name];
            if exist(fname,'file')
                [~,name,~] = fileparts(fname);
                name = strrep(name,'.nii','');
                break;
            end
        end
        if ~isempty(name), break; end
    end
    if isempty(name)
        disp(['input file is not found : ' argv]);
        return;
    end
    cachename = ['results/cache/n2r-info-' name '.mat'];
    if handles.noCache > 0 || ~exist(cachename,'file')
        info = niftiinfo(fname);
        if handles.noCache == 0
            save(cachename,'info');
        end
    else
        load(cachename);
    end

    % resampling volume
    step=size(atlasV,1)/info.ImageSize(1); stepZ=size(atlasV,3)/info.ImageSize(3);
    x2 = floor(size(atlasV,1)/step);
    y2 = floor(size(atlasV,2)/step);
    z2 = floor(size(atlasV,3)/stepZ);
    V2 = single(zeros(x2,y2,z2));
    for z=1:z2
        for y=1:y2
            for x=1:x2
                A = atlasV(x*step-(step-1):x*step,y*step-(step-1):y*step,z*stepZ-(stepZ-1):z*stepZ);
                m = mode(A(:));
                V2(x,y,z) = m;
            end
        end
    end
    atlasV = V2;

    % voxel indeces of 3D space
    idxs = cell(1,nodeNum);
    for i=1:nodeNum
        idx = find(atlasV==R(i));
        idxs{i} = idx;
    end

    % process each file
    CX = {}; CXm = {}; names = {}; savename = '';
    for i = 1:N
        % load multivariate time-series csv or mat file
        argv = handles.niiFiles{i};
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
            name = strrep(name,'.nii','');
            if isempty(savename)
                savename = name;
            end

            % read nifti file
            disp(['processing : ' name]);
            cachename = ['results/cache/n2r-' atlasName num2str(nodeNum) '-' name '.mat'];
            if handles.noCache > 0 || ~exist(cachename,'file')
                info = niftiinfo(fname);
                V = niftiread(info);
                V = adjustVolumeDir(V, info);

                if size(V,1) ~= size(atlasV,1) || size(V,2) ~= size(atlasV,2) || size(V,3) ~= size(atlasV,3)
                    disp(['atlas space and fMRI space does not match. ignoring : ' fname]);
                    continue;
                end

                X = single(zeros(nodeNum, size(V,4)));
                for t=1:size(V,4)
                    V2=squeeze(V(:,:,:,t));
                    for k=1:nodeNum
                        X(k,t) = nanmean(V2(idxs{k}));
                    end
                end
                if handles.noCache == 0
                    save(cachename,'X');
                end
            else
                load(cachename)
            end
            Xm = mean(X,2);
            X = X - Xm;

            % signal transform raw or not
            if handles.transform == 1
                [X, sig, c, maxsi, minsi] = convert2SigmoidSignal(X, handles.transopt);
            end

            % show output signals
            if handles.showSig > 0
                figure; plot(X.');
                title(['ROI Signals : ' strrep(name,'_','-')]);
                xlabel('Time Series');
                ylabel('Signal Value');
            end

            % show output signals
            if handles.showRas > 0
                figure; imagesc(X);
                title(['Raster plot of ROI Signals : ' strrep(name,'_','-')]);
                xlabel('Time Series');
                ylabel('Node number');
                colorbar;
            end

            CX{end+1} = X;
            CXm{end+1} = Xm;
            names{end+1} = name;
        end
    end
    
    if ~isempty(CX)
        saveResultFiles(handles, CX, CXm, names, savename);
    end
end

%%
% output result matrix files
%
%%
% output result matrix files
%
function saveResultFiles(handles, CX, CXm, names, outname)
    if handles.format == 1
        for i=1:length(CX)
            X = single(CX{i});
            m = CXm{i};
            outfname = [handles.outpath '/' names{i} '.mat'];
            save(outfname,'X','m');
            disp(['output mat file : ' outfname]);
        end
    elseif handles.format == 2
        outfname = [handles.outpath '/' outname '_all.mat'];
        save(outfname, 'CX', 'CXm', 'names', '-v7.3');
        disp(['output mat file : ' outfname]);
    else
        % output result matrix csv file
        for i=1:length(CX)
            X = single(CX{i});
            outputCsvFile(X, [handles.outpath '/' names{i} '.csv']);
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

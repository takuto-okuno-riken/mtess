%%
% !! full disk cache calculation to reduce memory consumption !!
% Caluclate MTESS (cross 2 groups), MTESS statistical properties, Node MTESS and Node MTESS statistical properties
% returns MTESS matrix (cell number x cell number)(MTS), MTESS statistical property matrix (cell number x cell number x 7)(MTSp), 
%   Node MTESS (cell number x cell number x node)(nMTS) and Node MTESS statistical properties (cell number x cell number x node x 7)(nMTSp).
%   Data in the middle of calculation, such as mean (Means), standard deviation (Stds), DFT amplitude (Amps), correlation matrix (FCs),
%   partial correlation matrix (PCs), cross-correlation matrix (CCs) and partial cross-correlation matrix (PCCs).
% input:
%  CX               cells of multivariate time series matrix {(node x time series)} x cell number (time series length can be different)
%  CY               cells of multivariate time series matrix {(node x time series)} x cell number (time series length can be different)
%  range            range [min, max] of time series for normalized mean and std dev (default: min and max of input CX)
%  nDFT             DFT sampling number (even number) (default: 100)
%  pccFunc          Partial Cross-Correlation function (default: @calcPartialCrossCorrelation)
%  ccLags           time lags for Cross-Correlation function (default: 8)
%  pccLags          time lags for Partial Cross-Correlation function (default: 8)
%  CXNames          CX signals names used for cache filename (default: {})
%  CYNames          CY signals names used for cache filename (default: {})
%  cachePath        cache file path (default: 'results/cache')

function [MTS, MTSp, nMTS, nMTSp] = calcMtessCross_c(CX, CY, range, nDft, pccFunc, ccLags, pccLags, CXNames, CYNames, cachePath)
    if nargin < 10, cachePath = 'results/cache'; end 
    if nargin < 9, CYNames = {}; end 
    if nargin < 8, CXNames = {}; end 
    if nargin < 7, pccLags = 8; end
    if nargin < 6, ccLags = 8; end
    if nargin < 5, pccFunc = @calcPartialCrossCorrelation; end
    if nargin < 4, nDft = 100; end
    if nargin < 3, range = NaN; end

    cxLen = length(CX);
    cyLen = length(CY);
    nodeNum = size(CX{1},1);
    memClass = class(CX{1});
    % check data file. node num should be same.
    for i=2:cxLen
        if size(CX{i},1) ~= nodeNum
            disp('Error : input time series should have same node number.');
            return;
        end
    end
    for i=1:cyLen
        if size(CY{i},1) ~= nodeNum
            disp('Error : input time series should have same node number.');
            return;
        end
    end

    % find time series range
    if isnan(range)
        minv = min(CX{1},[],'all');
        maxv = max(CX{1},[],'all');
        for i=2:cxLen
            v = min(CX{i},[],'all');
            if v < minv, minv = v; end
            v = max(CX{i},[],'all');
            if v > maxv, maxv = v; end
        end
        range = [minv, maxv];
    end
    tRange = range(2) - range(1); % should be positive

    % calc statistical properties
    if ~isempty(CXNames) && ~isempty(CYNames) && ~exist(cachePath,'dir')
        mkdir(cachePath);
    end
    if isequal(pccFunc,@calcSvPartialCrossCorrelation), palgo='svp';
    elseif isequal(pccFunc,@calcPcPartialCrossCorrelation), palgo='pcp';
    else, palgo='p';
    end

    calcStatProps(CX, nDft, pccFunc, palgo, ccLags, pccLags, CXNames, cachePath);
    calcStatProps(CY, nDft, pccFunc, palgo, ccLags, pccLags, CYNames, cachePath);

    % calc MTESS
    A = single(nan(nodeNum)); A= tril(A,0); % half does not support
    MTSp = nan(cxLen,cyLen,7,memClass);
    nMTSp = nan(cxLen,cyLen,nodeNum,7,memClass);
    for i=1:cxLen
        disp(['load cache of ' CXNames{i}]);
        [sz1, sz2] = size(CX{i});
        cachef = [cachePath '/mtess-' CXNames{i} '-' num2str(sz1) 'x' num2str(sz2) 'd' num2str(nDft) 'c' num2str(ccLags) palgo num2str(pccLags) '.mat'];
        fi = load(cachef);

        B = nan(cyLen,7,memClass);
        nB = nan(cyLen,nodeNum,7,memClass);
        parfor j=1:cyLen
            tj=tic;
            C = nan(7,1,memClass);
            nC = nan(nodeNum,7,memClass);
            [sz1, sz2] = size(CY{j});
            cachef = [cachePath '/mtess-' CYNames{j} '-' num2str(sz1) 'x' num2str(sz2) 'd' num2str(nDft) 'c' num2str(ccLags) palgo num2str(pccLags) '.mat'];
            fj = load(cachef);

            % calc mean distance (normalized)
            dm = fi.xm - fj.xm;
            nC(:,1) = abs(dm) / (tRange/2); % normalize

            % calc std deviation distance (normalized)
            ds = fi.xsd - fj.xsd;
            nC(:,2) = abs(ds) / (tRange/4); % normalize

            % calc amplitude similarity
            C(3) = 5 * getCosSimilarity(fi.xamp,fj.xamp);
            A1 = single(squeeze(fi.xamp));
            A2 = single(squeeze(fj.xamp));
            for k=1:nodeNum
                nC(k,3) = 5 * getCosSimilarity(A1(k,:),A2(k,:));
            end
            
            % calc zero-lag covariance similarity
            FC1 = single(squeeze(fi.xcc(:,:,ccLags+1)));
            FC2 = single(squeeze(fj.xcc(:,:,ccLags+1)));
            C(4) = 5 * getCosSimilarity(FC1+A, FC2+A); % half does not work
            for k=1:nodeNum
                nC(k,4) = 5 * getCosSimilarity([FC1(k,:)+A(k,:), (FC1(:,k)+A(:,k)).'], [FC2(k,:)+A(k,:), (FC2(:,k)+A(:,k)).']);
            end
            
            % calc zero-lag partial covariance similarity
            if ~isempty(fi.xpcc) && ~isempty(fj.xpcc)
                PC1 = single(squeeze(fi.xpcc(:,:,pccLags+1)));
                PC2 = single(squeeze(fj.xpcc(:,:,pccLags+1)));
            else
                PC1 = single(fi.pc);
                PC2 = single(fj.pc);
            end
            C(5) = 5 * getCosSimilarity(PC1+A, PC2+A); % half does not work
            for k=1:nodeNum
                nC(k,5) = 5 * getCosSimilarity([PC1(k,:)+A(k,:), (PC1(:,k)+A(:,k)).'], [PC2(k,:)+A(k,:), (PC2(:,k)+A(:,k)).']);
            end
            
            % calc cross-covariance simirality
            CC1 = single(squeeze(fi.xcc(:,:,[1:ccLags,ccLags+2:end])));
            CC2 = single(squeeze(fj.xcc(:,:,[1:ccLags,ccLags+2:end])));
            C(6) = 5 * getCosSimilarity(CC1+A,CC2+A); % half does not work
            for k=1:nodeNum
                R1 = [CC1(k,:,:)+A(k,:), permute(CC1(:,k,:)+A(:,k),[2 1 3])];
                R2 = [CC2(k,:,:)+A(k,:), permute(CC2(:,k,:)+A(:,k),[2 1 3])];
                nC(k,6) = 5 * getCosSimilarity(R1, R2);
            end

            % calc partial cross-covariance simirality
            if ~isempty(fi.xpcc) && ~isempty(fj.xpcc)
                PCC1 = single(squeeze(fi.xpcc(:,:,[1:pccLags,pccLags+2:end])));
                PCC2 = single(squeeze(fj.xpcc(:,:,[1:pccLags,pccLags+2:end])));
                C(7) = 5 * getCosSimilarity(PCC1+A,PCC2+A); % half does not work
                for k=1:nodeNum
                    R1 = [PCC1(k,:,:)+A(k,:), permute(PCC1(:,k,:)+A(:,k),[2 1 3])];
                    R2 = [PCC2(k,:,:)+A(k,:), permute(PCC2(:,k,:)+A(:,k),[2 1 3])];
                    nC(k,7) = 5 * getCosSimilarity(R1, R2);
                end
            end
            B(j,:) = C;
            nB(j,:,:) = nC;
            s = toc(tj);
            disp([num2str(i) '-' num2str(j) ' ' num2str(s) 'sec']);
        end
        MTSp(i,:,:) = B;
        nMTSp(i,:,:,:) = nB;
        clear fi;
        % memory leak in working pool?? shutdown it
%        poolobj = gcp('nocreate');
%        delete(poolobj);
    end

    % calc mean and std dev similarity
    m1 = 5 * (1 - nanmean(nMTSp(:,:,:,1),3));
    m2 = 5 * (1 - nanmean(nMTSp(:,:,:,2),3));
    m1(m1<0) = 0;
    m2(m2<0) = 0;
    MTSp(:,:,1) = m1;
    MTSp(:,:,2) = m2;
    m1 = 5 * (1 - nMTSp(:,:,:,1));
    m2 = 5 * (1 - nMTSp(:,:,:,2));
    m1(m1<0) = 0;
    m2(m2<0) = 0;
    nMTSp(:,:,:,1) = m1;
    nMTSp(:,:,:,2) = m2;

    % calc MTESS & Node MTESS
    MTSp(MTSp<0)=0;
    nMTSp(nMTSp<0)=0;
    MTSp(MTSp>5) = 5; % this may happen because of decimal point calculation
    nMTSp(nMTSp>5) = 5; % this may happen because of decimal point calculation
    MTS = nanmean(MTSp,3);
    nMTS = nanmean(nMTSp,4);
end

function calcStatProps(CX, nDft, pccFunc, palgo, ccLags, pccLags, CXNames, cachePath)
    cxLen = length(CX);
    for nn=1:cxLen
        X = CX{nn};
        if ~isempty(CXNames)
            cachef = [cachePath '/mtess-' CXNames{nn} '-' num2str(size(X,1)) 'x' num2str(size(X,2)) 'd' num2str(nDft) 'c' num2str(ccLags) palgo num2str(pccLags) '.mat'];
        end
        if ~isempty(CXNames) && exist(cachef,'file')
            disp(['cache exist : ' CXNames{nn}]);
        else
            xm = mean(X,2);
            xsd = std(X,1,2);
            xamp = calcDft(single(X),nDft); % half might take 'Inf'
            tc = tic;
            xcc = calcCrossCorrelation_(X,[],[],[],ccLags,0,false); % use gpu false
            s = toc(tc); disp([num2str(nn) ' : calcCrossCorr ' num2str(s) ' sec'])
            pc = [];
            tc = tic;
            if isempty(pccFunc)
                pc = calcPartialCorrelation_(X,[],[],[],0,false); % use gpu
                xpcc = [];
            elseif isequal(pccFunc,@calcSvPartialCrossCorrelation)
                xpcc = pccFunc(X,[],[],[],pccLags,'gaussian');
            else
                xpcc = pccFunc(X,[],[],[],pccLags,0,false); % use gpu
            end
            s = toc(tc); disp([num2str(nn) ' : calcPartialCrossCorr ' num2str(s) ' sec'])
            if ~isempty(CXNames)
                disp(['save cache of ' CXNames{nn}]);
                save(cachef, 'xm', 'xsd', 'xamp', 'xcc', 'xpcc', 'pc', '-v7.3');
            end
        end
    end
end

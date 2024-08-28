%%
% Caluclate MTESS, MTESS statistical properties, Node MTESS and Node MTESS statistical properties
% returns MTESS matrix (cell number x cell number)(MTS), MTESS statistical property matrix (cell number x cell number x 8)(MTSp), 
%   Node MTESS (cell number x cell number x node)(nMTS) and Node MTESS statistical properties (cell number x cell number x node x 8)(nMTSp).
%   Data in the middle of calculation, such as mean (Means), standard deviation (Stds), DFT amplitude (Amps), correlation matrix (FCs),
%   partial correlation matrix (PCs), cross-correlation matrix (CCs), partial cross-correlation matrix (PCCs) and multivariate kurtosis (mKTs).
% input:
%  CX               cells of multivariate time series matrix {(node x time series)} x cell number (time series length can be different)
%  range            range [min, max] of time series for normalized mean and std dev (default: min and max of input CX)
%  pccFunc          Partial Cross-Correlation function (default: @calcPartialCrossCorrelation)
%  acLags           time lags for Auto-Correlation function (default: 5)
%  pacLags          time lags for Partial Auto-Correlation function (default: 13)
%  ccLags           time lags for Cross-Correlation function (default: 2)
%  pccLags          time lags for Partial Cross-Correlation function (default: 4)
%  CXNames          CX signals names used for cache filename (default: {})
%  cachePath        cache file path (default: 'results/cache')

function [MTS, MTSp, nMTS, nMTSp, Means, Stds, ACs, PACs, FCs, PCs, CCs, PCCs, mKTs] = calcMtess(CX, range, pccFunc, acLags, pacLags, ccLags, pccLags, CXNames, cachePath)
    if nargin < 9, cachePath = 'results/cache'; end 
    if nargin < 8, CXNames = {}; end
    if nargin < 7, pccLags = 4; end
    if nargin < 6, ccLags = 2; end
    if nargin < 5, pacLags = 13; end
    if nargin < 4, acLags = 5; end
    if nargin < 3, pccFunc = @calcPartialCrossCorrelation; end
    if nargin < 2, range = NaN; end

    itemNum = 8;
    cLen = length(CX);
    nodeNum = size(CX{1},1);
    memClass = class(CX{1});
    % check data file. node num should be same.
    for i=2:cLen
        if size(CX{i},1) ~= nodeNum
            disp('Error : input time series should have same node number.');
            return;
        end
    end

    % find time series range
    if isnan(range)
        minv = min(CX{1},[],'all');
        maxv = max(CX{1},[],'all');
        for i=2:cLen
            v = min(CX{i},[],'all');
            if v < minv, minv = v; end
            v = max(CX{i},[],'all');
            if v > maxv, maxv = v; end
        end
        range = [minv, maxv];
    end
    tRange = range(2) - range(1); % should be positive

    % calc statistical properties
    if ~isempty(CXNames) && ~exist(cachePath,'dir')
        mkdir(cachePath);
    end
    if isequal(pccFunc,@calcSvPartialCrossCorrelation) palgo='svp';
    elseif isequal(pccFunc,@calcPcPartialCrossCorrelation) palgo='pcp';
    else palgo='p'; end

    Means = nan(cLen,nodeNum,memClass);
    Stds = nan(cLen,nodeNum,memClass);
    ACs = nan(cLen,nodeNum,acLags+1,memClass);
    PACs = nan(cLen,nodeNum,pacLags+1,memClass);
    FCs = nan(cLen,nodeNum,nodeNum,memClass);
    PCs = nan(cLen,nodeNum,nodeNum,memClass);
    CCs = nan(cLen,nodeNum,nodeNum,2*ccLags+1,memClass);
    if iscell(pccLags)
        lambda = pccLags{2}; % for redge regression
        pccLags = pccLags{1};
    else
        lambda = 0; % for redge regression
    end
    PCCs = nan(cLen,nodeNum,nodeNum,2*pccLags+1,memClass);
    mKTs = nan(cLen,1,memClass);
    for nn=1:cLen
        X = CX{nn};
        if ~isempty(CXNames)
            cachef = [cachePath '/mtess-' CXNames{nn} '-' num2str(size(X,1)) 'x' num2str(size(X,2)) 'a' num2str(acLags) 'c' num2str(ccLags) palgo num2str(pccLags) '.mat'];
        end
        if ~isempty(CXNames) && exist(cachef,'file')
            disp(['load cache of ' CXNames{nn}]);
            load(cachef);
        else
            xm = mean(X,2);
            xsd = std(X,1,2);
            xac = calcAutoCorrelation(X,acLags);
            xpac = calcPartialAutoCorrelation(X,pacLags);
            xcc = calcCrossCorrelation_(X,[],[],[],ccLags); % faster version
            if isempty(pccFunc)
                xpcc = [];
            elseif isequal(pccFunc,@calcSvPartialCrossCorrelation)
                xpcc = pccFunc(X,[],[],[],pccLags,'gaussian');
            elseif isequal(pccFunc,@calcPartialCrossCorrelation)
                xpcc = pccFunc(X,[],[],[],pccLags,0,lambda); % linear or ridge regress
            else
                xpcc = pccFunc(X,[],[],[],pccLags);
            end
            [~, xmkt] = calcMskewKurt(X);
            if ~isempty(CXNames)
                disp(['save cache of ' CXNames{nn}]);
                save(cachef, 'xm', 'xsd', 'xac', 'xpac', 'xcc', 'xpcc', 'xmkt');
            end
        end        
        Means(nn,:) = xm;
        Stds(nn,:) = xsd;
        ACs(nn,:,:) = xac;
        PACs(nn,:,:) = xpac;
        CCs(nn,:,:,:) = xcc;
        FCs(nn,:,:) = squeeze(xcc(:,:,ccLags+1));
        if ~isempty(xpcc)
            PCCs(nn,:,:,:) = xpcc;
            PCs(nn,:,:) = squeeze(xpcc(:,:,pccLags+1));
        end
        mKTs(nn) = xmkt;
    end

    % calc MTESS
    A = single(nan(nodeNum)); A= tril(A,0); % half does not support
    MTSp = nan(cLen,cLen,itemNum,memClass);
    nMTSp = nan(cLen,cLen,nodeNum,itemNum,memClass);
    for i=1:cLen
        B = nan(cLen,itemNum,memClass);
        nB = nan(cLen,nodeNum,itemNum,memClass);
        parfor j=i+1:cLen
            C = nan(itemNum,1,memClass);
            nC = nan(nodeNum,itemNum,memClass);

            % calc std deviation distance (normalized)
            ds = Stds(i,:)-Stds(j,:);
            D = 5 * (1 - abs(ds) / (tRange/4)); % normalize
            D(D<0) = 0;
            C(1) = nanmean(D);
            nC(:,1) = D;

            % calc AC/PAC similarity
            C(2) = 5 * getCosSimilarity(ACs(i,:,2:end),ACs(j,:,2:end));
            A1 = squeeze(ACs(i,:,:));
            A2 = squeeze(ACs(j,:,:));
            for k=1:nodeNum
                nC(k,2) = 5 * getCosSimilarity(A1(k,2:end),A2(k,2:end));
            end
            C(3) = 5 * getCosSimilarity(PACs(i,:,2:end),PACs(j,:,2:end));
            A1 = squeeze(PACs(i,:,:));
            A2 = squeeze(PACs(j,:,:));
            for k=1:nodeNum
                nC(k,3) = 5 * getCosSimilarity(A1(k,2:end),A2(k,2:end));
            end
            
            % calc zero-lag covariance similarity
            FC1 = squeeze(FCs(i,:,:)) + A;
            FC2 = squeeze(FCs(j,:,:)) + A;
            C(4) = 5 * getCosSimilarity(FC1, FC2);
            for k=1:nodeNum
                nC(k,4) = 5 * getCosSimilarity([FC1(k,:), FC1(:,k).'], [FC2(k,:), FC2(:,k).']);
            end
            
            % calc zero-lag partial covariance similarity
            if ~isempty(pccFunc)
                PC1 = squeeze(PCs(i,:,:)) + A;
                PC2 = squeeze(PCs(j,:,:)) + A;
                C(5) = 5 * getCosSimilarity(PC1, PC2);
                for k=1:nodeNum
                    nC(k,5) = 5 * getCosSimilarity([PC1(k,:), PC1(:,k).'], [PC2(k,:), PC2(:,k).']);
                end
            end
            
            % calc cross-covariance simirality
            CC1 = squeeze(CCs(i,:,:,[1:ccLags,ccLags+2:end])) + A;
            CC2 = squeeze(CCs(j,:,:,[1:ccLags,ccLags+2:end])) + A;
            C(6) = 5 * getCosSimilarity(CC1,CC2);
            for k=1:nodeNum
                R1 = [CC1(k,:,:), permute(CC1(:,k,:),[2 1 3])];
                R2 = [CC2(k,:,:), permute(CC2(:,k,:),[2 1 3])];
                nC(k,6) = 5 * getCosSimilarity(R1, R2);
            end

            % calc partial cross-covariance simirality
            if ~isempty(pccFunc)
                PCC1 = squeeze(PCCs(i,:,:,[1:pccLags,pccLags+2:end])) + A;
                PCC2 = squeeze(PCCs(j,:,:,[1:pccLags,pccLags+2:end])) + A;
                C(7) = 5 * getCosSimilarity(PCC1,PCC2);
                for k=1:nodeNum
                    R1 = [PCC1(k,:,:), permute(PCC1(:,k,:),[2 1 3])];
                    R2 = [PCC2(k,:,:), permute(PCC2(:,k,:),[2 1 3])];
                    nC(k,7) = 5 * getCosSimilarity(R1, R2);
                end
            end

            % multivariate kurtosis
            mkt = nodeNum*(nodeNum+2) / 2; % 2 is empirically defined.
            ds = 5 * (1 - abs(mKTs(i)-mKTs(j)) / mkt); % normalize
            if ds < 0, ds = 0; end
            C(8) = ds;
            nC(:,8) = ds;
            % output
            B(j,:) = C;
            nB(j,:,:) = nC;
        end
        MTSp(i,:,:) = B;
        nMTSp(i,:,:,:) = nB;
    end

    % calc MTESS & Node MTESS
    MTSp(MTSp<0)=0;
    nMTSp(nMTSp<0)=0;
    MTSp(MTSp>5) = 5; % this may happen because of decimal point calculation
    nMTSp(nMTSp>5) = 5; % this may happen because of decimal point calculation
    MTS = nanmean(MTSp,3);
    nMTS = nanmean(nMTSp,4);
end

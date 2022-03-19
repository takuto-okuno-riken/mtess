%%
% Caluclate MTESS, MTESS statistical properties, Node MTESS and Node MTESS statistical properties
% returns MTESS matrix (cell number x cell number)(MTS), MTESS statistical property matrix (cell number x cell number x 7)(MTSp), 
%   Node MTESS (cell number x cell number x node)(nMTS) and Node MTESS statistical properties (cell number x cell number x node x 7)(nMTSp).
%   Data in the middle of calculation, such as mean (Means), standard deviation (Stds), DFT amplitude (Amps), correlation matrix (FCs),
%   partial correlation matrix (PCs), cross-correlation matrix (CCs) and partial cross-correlation matrix (PCCs).
% input:
%  CX               cells of multivariate time series matrix {(node x time series)} x cell number (time series length can be different)
%  range            range [min, max] of time series for normalized mean and std dev (default: min and max of input CX)
%  nDFT             DFT sampling number (even number) (default: 100)
%  pccFunc          Partial Cross-Correlation function (default: @calcPartialCrossCorrelation)
%  ccLags           time lags for Cross-Correlation function (default: 8)
%  pccLags          time lags for Partial Cross-Correlation function (default: 8)
%  CXNames          CX signals names used for cache filename (default: {})

function [MTS, MTSp, nMTS, nMTSp, Means, Stds, Amps, FCs, PCs, CCs, PCCs] = calcMtess(CX, range, nDft, pccFunc, ccLags, pccLags, CXNames)
    if nargin < 7, CXNames = {}; end
    if nargin < 6, pccLags = 8; end
    if nargin < 5, ccLags = 8; end
    if nargin < 4, pccFunc = @calcPartialCrossCorrelation; end
    if nargin < 3, nDft = 100; end
    if nargin < 2, range = NaN; end

    cLen = length(CX);
    nodeNum = size(CX{1},1);
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
    if ~isempty(CXNames) && ~exist('results/cache','dir')
        mkdir('results/cache');
    end
    if isequal(pccFunc,@calcSvPartialCrossCorrelation) palgo='svp';
    elseif isequal(pccFunc,@calcPcPartialCrossCorrelation) palgo='pcp';
    else palgo='p'; end

    Means = single(nan(cLen,nodeNum));
    Stds = single(nan(cLen,nodeNum));
    Amps = single(nan(cLen,nodeNum,nDft/2-1));
    FCs = single(nan(cLen,nodeNum,nodeNum));
    PCs = single(nan(cLen,nodeNum,nodeNum));
    CCs = single(nan(cLen,nodeNum,nodeNum,2*ccLags+1));
    PCCs = single(nan(cLen,nodeNum,nodeNum,2*pccLags+1));
    for nn=1:cLen
        X = CX{nn};
        if ~isempty(CXNames)
            cachef = ['results/cache/mtess-' CXNames{nn} '-' num2str(size(X,1)) 'x' num2str(size(X,2)) 'd' num2str(nDft) 'c' num2str(ccLags) palgo num2str(pccLags) '.mat'];
        end
        if ~isempty(CXNames) && exist(cachef,'file')
            disp(['load cache of ' CXNames{nn}]);
            load(cachef);
        else
            xm = single(mean(X,2));
            xsd = single(std(X,1,2));
            xamp = single(calcDft(X,nDft));
            xcc = single(calcCrossCorrelation(X,[],[],[],ccLags));
            if isequal(pccFunc,@calcSvPartialCrossCorrelation)
                xpcc = single(pccFunc(X,[],[],[],pccLags,'gaussian'));
            else
                xpcc = single(pccFunc(X,[],[],[],pccLags));
            end
            if ~isempty(CXNames)
                disp(['save cache of ' CXNames{nn}]);
                save(cachef, 'xm', 'xsd', 'xamp', 'xcc', 'xpcc');
            end
        end        
        Means(nn,:) = xm;
        Stds(nn,:) = xsd;
        Amps(nn,:,:) = xamp;
        CCs(nn,:,:,:) = xcc;
        PCCs(nn,:,:,:) = xpcc;
        FCs(nn,:,:) = squeeze(xcc(:,:,ccLags+1));
        PCs(nn,:,:) = squeeze(xpcc(:,:,pccLags+1));
%        PCs(nn,:,:) = calcPartialCorrelation(si);
    end
    Amps = Amps * (nDft/2) / (tRange/2); % normalize

    % calc MTESS
    A = single(nan(nodeNum)); A= tril(A,0);
    MTSp = single(nan(cLen,cLen,7));
    nMTSp = single(nan(cLen,cLen,nodeNum,7));
    for i=1:cLen
        for j=i+1:cLen
            % calc mean distance (normalized)
            dm = Means(i,:)-Means(j,:);
            nMTSp(i,j,:,1) = abs(dm) / (tRange/2); % normalize

            % calc std deviation distance (normalized)
            ds = Stds(i,:)-Stds(j,:);
            nMTSp(i,j,:,2) = abs(ds) / (tRange/4); % normalize

            % calc amplitude similarity
            MTSp(i,j,3) = 5 * getCosSimilarity(Amps(i,:,:),Amps(j,:,:));
            for k=1:nodeNum
                nMTSp(i,j,k,3) = 5 * getCosSimilarity(Amps(i,k,:),Amps(j,k,:));
            end
            
            % calc zero-lag covariance similarity
            FC1 = squeeze(FCs(i,:,:));
            FC2 = squeeze(FCs(j,:,:));
            MTSp(i,j,4) = 5 * getCosSimilarity(FC1+A, FC2+A);
            for k=1:nodeNum
                nMTSp(i,j,k,4) = 5 * getCosSimilarity([FC1(k,:)+A(k,:), (FC1(:,k)+A(:,k)).'], [FC2(k,:)+A(k,:), (FC2(:,k)+A(:,k)).']);
            end
            
            % calc zero-lag partial covariance similarity
            PC1 = squeeze(PCs(i,:,:));
            PC2 = squeeze(PCs(j,:,:));
            MTSp(i,j,5) = 5 * getCosSimilarity(PC1+A, PC2+A);
            for k=1:nodeNum
                nMTSp(i,j,k,5) = 5 * getCosSimilarity([PC1(k,:)+A(k,:), (PC1(:,k)+A(:,k)).'], [PC2(k,:)+A(k,:), (PC2(:,k)+A(:,k)).']);
            end
            
            % calc cross-covariance simirality
            CC1 = squeeze(CCs(i,:,:,[1:ccLags,ccLags+2:end]));
            CC2 = squeeze(CCs(j,:,:,[1:ccLags,ccLags+2:end]));
            MTSp(i,j,6) = 5 * getCosSimilarity(CC1+A,CC2+A);
            for k=1:nodeNum
                R1 = [CC1(k,:,:)+A(k,:), permute(CC1(:,k,:)+A(:,k),[2 1 3])];
                R2 = [CC2(k,:,:)+A(k,:), permute(CC2(:,k,:)+A(:,k),[2 1 3])];
                nMTSp(i,j,k,6) = 5 * getCosSimilarity(R1, R2);
            end

            % calc partial cross-covariance simirality
            PCC1 = squeeze(PCCs(i,:,:,[1:pccLags,pccLags+2:end]));
            PCC2 = squeeze(PCCs(j,:,:,[1:pccLags,pccLags+2:end]));
            MTSp(i,j,7) = 5 * getCosSimilarity(PCC1+A,PCC2+A);
            for k=1:nodeNum
                R1 = [PCC1(k,:,:)+A(k,:), permute(PCC1(:,k,:)+A(:,k),[2 1 3])];
                R2 = [PCC2(k,:,:)+A(k,:), permute(PCC2(:,k,:)+A(:,k),[2 1 3])];
                nMTSp(i,j,k,7) = 5 * getCosSimilarity(R1, R2);
            end
        end
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

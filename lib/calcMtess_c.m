%%
% !! full disk cache calculation to reduce memory consumption !!
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
%  cachePath        cache file path (default: 'results/cache')
%  pccPca           PCA structure for pcc calculation (default: [])

function [MTS, MTSp, nMTS, nMTSp] = calcMtess_c(CX, range, nDft, pccFunc, ccLags, pccLags, CXNames, cachePath, pccPca)
    if nargin < 9, pccPca = []; end
    if nargin < 8, cachePath = 'results/cache'; end 
    if nargin < 7, CXNames = {}; end 
    if nargin < 6, pccLags = 8; end
    if nargin < 5, ccLags = 8; end
    if nargin < 4, pccFunc = @calcPartialCrossCorrelation; end
    if nargin < 3, nDft = 100; end
    if nargin < 2, range = NaN; end

    cLen = length(CX);
    nodeNum = size(CX{1},1);
    memClass = class(CX{1});
    isnMTS = nargout >= 3;
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
    if isequal(pccFunc,@calcSvPartialCrossCorrelation), palgo='svp';
    elseif isequal(pccFunc,@calcPcPartialCrossCorrelation), palgo='pcp';
    else, palgo='p';
    end
    if isnMTS, ostr=''; else, ostr='n'; end
    A = ones(nodeNum,'logical'); A = triu(A,1);
    if isempty(pccPca)
        B = A;
    else
        B = ones(pccPca.maxComp,'logical'); B = triu(B,1);
    end

    for nn=1:cLen
        X = CX{nn};
        if ~isempty(CXNames)
            cachef = [cachePath '/mtess-' CXNames{nn} '-' num2str(size(X,1)) 'x' num2str(size(X,2)) 'd' num2str(nDft) 'c' num2str(ccLags) palgo num2str(pccLags) ostr '.mat'];
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
            if ~isnMTS && ~isempty(xcc)  % no node MTESS output mode
                lsz = size(xcc,3);
                xcc2 = zeros(nodeNum*(nodeNum-1)/2,lsz,memClass);
                for i=1:lsz
                    T = xcc(:,:,i);
                    xcc2(:,i) = T(A);
                end
                xcc = xcc2; clear xcc2;
            end
            pc = [];
            tc = tic;
            if ~isempty(pccPca)
                X = (X' - pccPca.mu) * pccPca.invCoeff;
                X = X.';
            end
            if isempty(pccFunc)
                xpcc = [];
            elseif isequal(pccFunc,@calcPartialCorrelation) || isequal(pccFunc,@calcPartialCorrelation_) || isequal(pccFunc,@calcPartialCorrelation__)
                pc = pccFunc(X,[],[],[],0);
            elseif isequal(pccFunc,@calcSvPartialCrossCorrelation)
                xpcc = pccFunc(X,[],[],[],pccLags,'gaussian');
            else
                xpcc = pccFunc(X,[],[],[],pccLags,0,false); % use gpu
            end
            s = toc(tc); disp([num2str(nn) ' : calcPartialCrossCorr ' num2str(s) ' sec'])
            if ~isnMTS && ~isempty(xpcc) % no node MTESS output mode
                [sz1, ~, sz3] = size(xpcc);
                xpcc2 = zeros(sz1*(sz1-1)/2,sz3,memClass);
                for i=1:sz3
                    T = xpcc(:,:,i);
                    xpcc2(:,i) = T(B);
                end
                xpcc = xpcc2; clear xpcc2;
            end
            if ~isempty(CXNames)
                disp(['save cache of ' CXNames{nn}]);
                save(cachef, 'xm', 'xsd', 'xamp', 'xcc', 'xpcc', 'pc', '-v7.3');
            end
        end
    end

    % calc MTESS
    A = nan(nodeNum,'single'); A= tril(A,0); % half does not support
    [sz1, sz2] = size(CX{1});
    procf = [cachePath '/mtess-' CXNames{1} '-' num2str(sz1) 'x' num2str(sz2) 'd' num2str(nDft) 'c' num2str(ccLags) palgo num2str(pccLags) ostr '-proc.mat'];
    if ~isempty(CXNames) && exist(procf,'file')
        load(procf);
    else
        MTSp = nan(cLen,cLen,7,memClass);
        if isnMTS
            nMTSp = nan(cLen,cLen,nodeNum,7,memClass);
        else
            nMTSp = [];
        end
    end
    for i=1:cLen
        if ~isnan(MTSp(i,cLen,3))
            disp(['skip : ' CXNames{i}]);
            continue;
        end
        disp(['load cache of ' CXNames{i}]);
        [sz1, sz2] = size(CX{i});
        cachef = [cachePath '/mtess-' CXNames{i} '-' num2str(sz1) 'x' num2str(sz2) 'd' num2str(nDft) 'c' num2str(ccLags) palgo num2str(pccLags) ostr '.mat'];
        fi = load(cachef);

        B = nan(cLen,7,memClass);
        if isnMTS, nB = nan(cLen,nodeNum,7,memClass); end
        parfor j=i+1:cLen
            tj=tic;
            C = nan(7,1,memClass);
            if isnMTS, nC = nan(nodeNum,7,memClass); end
            [sz1, sz2] = size(CX{j});
            cachef = [cachePath '/mtess-' CXNames{j} '-' num2str(sz1) 'x' num2str(sz2) 'd' num2str(nDft) 'c' num2str(ccLags) palgo num2str(pccLags) ostr '.mat'];
            fj = load(cachef);

            % calc mean distance (normalized)
            dm = fi.xm - fj.xm;
            D = abs(dm) / (tRange/2); % normalize
            D(D<0) = 0;
            C(1) = 5 * (1 - nanmean(D));
            if isnMTS, nC(:,1)=D; end

            % calc std deviation distance (normalized)
            ds = fi.xsd - fj.xsd;
            D = abs(ds) / (tRange/4); % normalize
            D(D<0) = 0;
            C(2) = 5 * (1 - nanmean(D));
            if isnMTS, nC(:,2)=D; end

            % calc amplitude similarity
            C(3) = 5 * getCosSimilarity(fi.xamp,fj.xamp);
            if isnMTS
                A1 = single(squeeze(fi.xamp));
                A2 = single(squeeze(fj.xamp));
                for k=1:nodeNum
                    nC(k,3) = 5 * getCosSimilarity(A1(k,:),A2(k,:));
                end
            end            
            % calc zero-lag covariance similarity
            if isnMTS
                FC1 = single(squeeze(fi.xcc(:,:,ccLags+1)));
                FC2 = single(squeeze(fj.xcc(:,:,ccLags+1)));
                C(4) = 5 * getCosSimilarity(FC1+A, FC2+A); % half does not work
                for k=1:nodeNum
                    nC(k,4) = 5 * getCosSimilarity([FC1(k,:)+A(k,:), (FC1(:,k)+A(:,k)).'], [FC2(k,:)+A(k,:), (FC2(:,k)+A(:,k)).']);
                end
            else
                FC1 = single(squeeze(fi.xcc(:,ccLags+1)));
                FC2 = single(squeeze(fj.xcc(:,ccLags+1)));
                C(4) = 5 * getCosSimilarity(FC1, FC2); % half does not work
            end
            % calc zero-lag partial covariance similarity
            if ~isempty(fi.xpcc)
                if isnMTS
                    PC1 = single(squeeze(fi.xpcc(:,:,pccLags+1)));
                else
                    PC1 = single(squeeze(fi.xpcc(:,pccLags+1)));
                end
            elseif ~isempty(fi.pc)
                PC1 = single(fi.pc);
            else
                PC1 = [];
            end
            if ~isempty(fj.xpcc)
                if isnMTS
                    PC2 = single(squeeze(fj.xpcc(:,:,pccLags+1)));
                else
                    PC2 = single(squeeze(fj.xpcc(:,pccLags+1)));
                end
            elseif ~isempty(fj.pc)
                PC2 = single(fj.pc);
            else
                PC2 = [];
            end
            if ~isempty(PC1) && ~isempty(PC2)
                if isnMTS
                    C(5) = 5 * getCosSimilarity(PC1+A, PC2+A); % half does not work
                    for k=1:nodeNum
                        nC(k,5) = 5 * getCosSimilarity([PC1(k,:)+A(k,:), (PC1(:,k)+A(:,k)).'], [PC2(k,:)+A(k,:), (PC2(:,k)+A(:,k)).']);
                    end
                else
                    C(5) = 5 * getCosSimilarity(PC1, PC2); % half does not work
                end
            end

            % calc cross-covariance simirality
            if isnMTS
                CC1 = single(squeeze(fi.xcc(:,:,[1:ccLags,ccLags+2:end])));
                CC2 = single(squeeze(fj.xcc(:,:,[1:ccLags,ccLags+2:end])));
                C(6) = 5 * getCosSimilarity(CC1+A,CC2+A); % half does not work
                for k=1:nodeNum
                    R1 = [CC1(k,:,:)+A(k,:), permute(CC1(:,k,:)+A(:,k),[2 1 3])];
                    R2 = [CC2(k,:,:)+A(k,:), permute(CC2(:,k,:)+A(:,k),[2 1 3])];
                    nC(k,6) = 5 * getCosSimilarity(R1, R2);
                end
            else
                CC1 = single(squeeze(fi.xcc(:,[1:ccLags,ccLags+2:end])));
                CC2 = single(squeeze(fj.xcc(:,[1:ccLags,ccLags+2:end])));
                C(6) = 5 * getCosSimilarity(CC1,CC2); % half does not work
            end
            % calc partial cross-covariance simirality
            if ~isempty(fi.xpcc) && ~isempty(fj.xpcc)
                if isnMTS
                    PCC1 = single(squeeze(fi.xpcc(:,:,[1:pccLags,pccLags+2:end])));
                    PCC2 = single(squeeze(fj.xpcc(:,:,[1:pccLags,pccLags+2:end])));
                    C(7) = 5 * getCosSimilarity(PCC1+A,PCC2+A); % half does not work
                    for k=1:nodeNum
                        R1 = [PCC1(k,:,:)+A(k,:), permute(PCC1(:,k,:)+A(:,k),[2 1 3])];
                        R2 = [PCC2(k,:,:)+A(k,:), permute(PCC2(:,k,:)+A(:,k),[2 1 3])];
                        nC(k,7) = 5 * getCosSimilarity(R1, R2);
                    end
                else
                    PCC1 = single(squeeze(fi.xpcc(:,[1:pccLags,pccLags+2:end])));
                    PCC2 = single(squeeze(fj.xpcc(:,[1:pccLags,pccLags+2:end])));
                    C(7) = 5 * getCosSimilarity(PCC1,PCC2); % half does not work
                end
            end
            B(j,:) = C;
            if isnMTS, nB(j,:,:) = nC; end
            s = toc(tj);
            disp([num2str(i) '-' num2str(j) ' ' num2str(s) 'sec']);
        end
        MTSp(i,:,:) = B;
        if isnMTS, nMTSp(i,:,:,:) = nB; end
        clear fi;
        if ~isempty(CXNames) && (sz1 > 1500 || cLen > 100)   % if node number is large
            save(procf,'MTSp','nMTSp','-v7.3');
        end
        % memory leak in working pool?? shutdown it
        if isnMTS
%            poolobj = gcp('nocreate');
%            delete(poolobj);
        end
    end

    % calc mean and std dev similarity
    if isnMTS
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
    end
    % calc MTESS & Node MTESS
    MTSp(MTSp<0)=0;
    nMTSp(nMTSp<0)=0;
    MTSp(MTSp>5) = 5; % this may happen because of decimal point calculation
    nMTSp(nMTSp>5) = 5; % this may happen because of decimal point calculation
    MTS = nanmean(MTSp,3);
    nMTS = nanmean(nMTSp,4);
end

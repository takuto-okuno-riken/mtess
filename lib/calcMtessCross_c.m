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
%  pccFunc          Partial Cross-Correlation function (default: @calcPartialCrossCorrelation)
%  acLags           time lags for Auto-Correlation / Partial Auto-Correlation function (default: 15)
%  ccLags           time lags for Cross-Correlation function (default: 8)
%  pccLags          time lags for Partial Cross-Correlation function (default: 8)
%  CXNames          CX signals names used for cache filename (default: {})
%  CYNames          CY signals names used for cache filename (default: {})
%  cachePath        cache file path (default: 'results/cache')

function [MTS, MTSp, nMTS, nMTSp] = calcMtessCross_c(CX, CY, range, pccFunc, acLags, ccLags, pccLags, CXNames, CYNames, cachePath)
    if nargin < 10, cachePath = 'results/cache'; end 
    if nargin < 9, CYNames = {}; end 
    if nargin < 8, CXNames = {}; end 
    if nargin < 7, pccLags = 8; end
    if nargin < 6, ccLags = 8; end
    if nargin < 5, acLags = 15; end
    if nargin < 4, pccFunc = @calcPartialCrossCorrelation; end
    if nargin < 3, range = NaN; end

    itemNum = 8;
    cxLen = length(CX);
    cyLen = length(CY);
    nodeNum = size(CX{1},1);
    memClass = class(CX{1});
    isnMTS = nargout >= 3;
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
    if isnMTS, ostr=''; else, ostr='n'; end

    calcStatProps(CX, pccFunc, palgo, acLags, ccLags, pccLags, CXNames, cachePath, isnMTS, ostr);
    calcStatProps(CY, pccFunc, palgo, acLags, ccLags, pccLags, CYNames, cachePath, isnMTS, ostr);

    % calc MTESS
    A = nan(nodeNum,'single'); A= tril(A,0); % half does not support
    [sz1, sz2] = size(CX{1});
    procf = [cachePath '/mtess-' CXNames{1} '-' CYNames{1} '-' num2str(sz1) 'x' num2str(sz2) 'a' num2str(acLags) 'c' num2str(ccLags) palgo num2str(pccLags) ostr '-xproc.mat'];
    if ~isempty(CXNames) && exist(procf,'file')
        load(procf);
    else
        MTSp = nan(cxLen,cyLen,itemNum,memClass);
        if isnMTS
            nMTSp = nan(cxLen,cyLen,nodeNum,itemNum,memClass);
        else
            nMTSp = [];
        end
    end
    for i=1:cxLen
        if ~isnan(MTSp(i,cyLen,3))
            disp(['skip : ' CXNames{i}]);
            continue;
        end
        disp(['load cache of ' CXNames{i}]);
        [sz1, sz2] = size(CX{i});
        cachef = [cachePath '/mtess-' CXNames{i} '-' num2str(sz1) 'x' num2str(sz2) 'a' num2str(acLags) 'c' num2str(ccLags) palgo num2str(pccLags) ostr '.mat'];
        fi = load(cachef);

        B = nan(cyLen,itemNum,memClass);
        if isnMTS, nB = nan(cyLen,nodeNum,itemNum,memClass); end
        parfor j=1:cyLen
            tj=tic;
            C = nan(itemNum,1,memClass);
            if isnMTS, nC = nan(nodeNum,itemNum,memClass); end
            [sz1, sz2] = size(CY{j});
            cachef = [cachePath '/mtess-' CYNames{j} '-' num2str(sz1) 'x' num2str(sz2) 'a' num2str(acLags) 'c' num2str(ccLags) palgo num2str(pccLags) ostr '.mat'];
            fj = load(cachef);

            % calc std deviation distance (normalized)
            ds = fi.xsd - fj.xsd;
            D = 5 * (1 - abs(ds) / (tRange/4)); % normalize
            D(D<0) = 0;
            C(1) = nanmean(D);
            if isnMTS, nC(:,1)=D; end

            % calc AC/PAC similarity
            C(2) = 5 * getCosSimilarity(fi.xac,fj.xac);
            if isnMTS
                A1 = single(squeeze(fi.xac));
                A2 = single(squeeze(fj.xac));
                for k=1:nodeNum
                    nC(k,2) = 5 * getCosSimilarity(A1(k,2:end),A2(k,2:end));
                end
            end
            C(3) = 5 * getCosSimilarity(fi.xpac,fj.xpac);
            if isnMTS
                A1 = squeeze(fi.xpac);
                A2 = squeeze(fj.xpac);
                for k=1:nodeNum
                    nC(k,3) = 5 * getCosSimilarity(A1(k,2:end),A2(k,2:end));
                end
            end
            % calc zero-lag covariance similarity
            if isnMTS
                FC1 = single(squeeze(fi.xcc(:,:,ccLags+1)) + A);
                FC2 = single(squeeze(fj.xcc(:,:,ccLags+1)) + A);
                C(4) = 5 * getCosSimilarity(FC1, FC2); % half does not work
                for k=1:nodeNum
                    nC(k,4) = 5 * getCosSimilarity([FC1(k,:), FC1(:,k).'], [FC2(k,:), FC2(:,k).']);
                end
            else
                FC1 = single(squeeze(fi.xcc(:,ccLags+1)) + A);
                FC2 = single(squeeze(fj.xcc(:,ccLags+1)) + A);
                C(4) = 5 * getCosSimilarity(FC1, FC2); % half does not work
            end
            % calc zero-lag partial covariance similarity
            if ~isempty(fi.xpcc)
                if isnMTS
                    PC1 = single(squeeze(fi.xpcc(:,:,pccLags+1)) + A);
                else
                    PC1 = single(squeeze(fi.xpcc(:,pccLags+1)) + A);
                end
            elseif ~isempty(fi.pc)
                PC1 = single(fi.pc);
            else
                PC1 = [];
            end
            if ~isempty(fj.xpcc)
                if isnMTS
                    PC2 = single(squeeze(fj.xpcc(:,:,pccLags+1)) + A);
                else
                    PC2 = single(squeeze(fj.xpcc(:,pccLags+1)) + A);
                end
            elseif ~isempty(fj.pc)
                PC2 = single(fj.pc);
            else
                PC2 = [];
            end
            if ~isempty(PC1) && ~isempty(PC2)
                if isnMTS
                    C(5) = 5 * getCosSimilarity(PC1, PC2); % half does not work
                    for k=1:nodeNum
                        nC(k,5) = 5 * getCosSimilarity([PC1(k,:), PC1(:,k).'], [PC2(k,:), PC2(:,k).']);
                    end
                else
                    C(5) = 5 * getCosSimilarity(PC1, PC2); % half does not work
                end
            end

            % calc cross-covariance simirality
            if isnMTS
                CC1 = single(squeeze(fi.xcc(:,:,[1:ccLags,ccLags+2:end])) + A);
                CC2 = single(squeeze(fj.xcc(:,:,[1:ccLags,ccLags+2:end])) + A);
                C(6) = 5 * getCosSimilarity(CC1+A,CC2+A); % half does not work
                for k=1:nodeNum
                    R1 = [CC1(k,:,:), permute(CC1(:,k,:),[2 1 3])];
                    R2 = [CC2(k,:,:), permute(CC2(:,k,:),[2 1 3])];
                    nC(k,6) = 5 * getCosSimilarity(R1, R2);
                end
            else
                CC1 = single(squeeze(fi.xcc(:,[1:ccLags,ccLags+2:end])) + A);
                CC2 = single(squeeze(fj.xcc(:,[1:ccLags,ccLags+2:end])) + A);
                C(6) = 5 * getCosSimilarity(CC1,CC2); % half does not work
            end

            % calc partial cross-covariance simirality
            if ~isempty(fi.xpcc) && ~isempty(fj.xpcc)
                if isnMTS
                    PCC1 = single(squeeze(fi.xpcc(:,:,[1:pccLags,pccLags+2:end])) + A);
                    PCC2 = single(squeeze(fj.xpcc(:,:,[1:pccLags,pccLags+2:end])) + A);
                    C(7) = 5 * getCosSimilarity(PCC1,PCC2); % half does not work
                    for k=1:nodeNum
                        R1 = [PCC1(k,:,:), permute(PCC1(:,k,:),[2 1 3])];
                        R2 = [PCC2(k,:,:), permute(PCC2(:,k,:),[2 1 3])];
                        nC(k,7) = 5 * getCosSimilarity(R1, R2);
                    end
                else
                    PCC1 = single(squeeze(fi.xpcc(:,[1:pccLags,pccLags+2:end])) + A);
                    PCC2 = single(squeeze(fj.xpcc(:,[1:pccLags,pccLags+2:end])) + A);
                    C(7) = 5 * getCosSimilarity(PCC1,PCC2); % half does not work
                end
            end
            % multivariate kurtosis
            mkt = nodeNum*(nodeNum+2) / 2; % 2 is empirically defined.
            ds = 5 * (1 - abs(fi.xmkt - fj.xmkt) / mkt); % normalize
            if ds < 0, ds = 0; end
            C(8) = ds;
            if isnMTS, nC(:,8) = ds; end

            B(j,:) = C;
            if isnMTS, nB(j,:,:) = nC; end
            s = toc(tj);
            disp([num2str(i) '-' num2str(j) ' ' num2str(s) 'sec']);
        end
        MTSp(i,:,:) = B;
        if isnMTS, nMTSp(i,:,:,:) = nB; end
        clear fi;
        if ~isempty(CXNames) && (sz1 > 1500 || cxLen > 100 || cyLen > 100)   % if node number is large
            save(procf,'MTSp','nMTSp','-v7.3');
        end
        % memory leak in working pool?? shutdown it
        if isnMTS
%            poolobj = gcp('nocreate');
%            delete(poolobj);
        end
    end

    % calc MTESS & Node MTESS
    MTSp(MTSp<0)=0;
    nMTSp(nMTSp<0)=0;
    MTSp(MTSp>5) = 5; % this may happen because of decimal point calculation
    nMTSp(nMTSp>5) = 5; % this may happen because of decimal point calculation
    MTS = nanmean(MTSp,3);
    nMTS = nanmean(nMTSp,4);
end

function calcStatProps(CX, pccFunc, palgo, acLags, ccLags, pccLags, CXNames, cachePath, isnMTS, ostr)
    nodeNum = size(CX{1},1);
    memClass = class(CX{1});
    A = ones(nodeNum,'logical'); A = triu(A,1);
    cxLen = length(CX);
    for nn=1:cxLen
        X = CX{nn};
        if ~isempty(CXNames)
            cachef = [cachePath '/mtess-' CXNames{nn} '-' num2str(size(X,1)) 'x' num2str(size(X,2)) 'a' num2str(acLags) 'c' num2str(ccLags) palgo num2str(pccLags) ostr '.mat'];
        end
        if ~isempty(CXNames) && exist(cachef,'file')
            disp(['cache exist : ' CXNames{nn}]);
        else
            xm = mean(X,2);
            xsd = std(X,1,2);
            xac = calcAutoCorrelation(X,acLags);
            xpac = calcPartialAutoCorrelation(X,acLags);
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
                    xpcc2(:,i) = T(A);
                end
                xpcc = xpcc2; clear xpcc2;
            end
            [~, xmkt] = calcMskewKurt(X);
            if ~isempty(CXNames)
                disp(['save cache of ' CXNames{nn}]);
                save(cachef, 'xm', 'xsd', 'xac', 'xpac', 'xcc', 'xpcc', 'xmkt', 'pc', '-v7.3');
            end
        end
    end
end

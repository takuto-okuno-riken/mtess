
function testSurrogateMVAR
    % load signals
    load('data/ww32-1.mat');

    %% test pattern 1 
    for lags=1:5
        netMVAR = initMvarNetwork(X, [], [], [], lags);
        Y1 = surrogateMVAR(X, [], [], [], netMVAR, 'gaussian',10);
        Y2 = surrogateMVAR(X, [], [], [], netMVAR, 'residuals',10);
        Y3 = surrogateMFT(X,10);

        figure; plotTwoSignals(X,Y1(:,:,3));
        figure; plotTwoSignals(X,Y2(:,:,3));
        figure; plotTwoSignals(X,Y3(:,:,3));

        if Y1(:,:,3)==Y1(:,:,2), disp('bad surrogate : Y1'); end
        if Y2(:,:,3)==Y2(:,:,2), disp('bad surrogate : Y2'); end
        if Y3(:,:,3)==Y3(:,:,2), disp('bad surrogate : Y3'); end

        % similarity check
        S = X;
        S(:,:,2) = X;
        checkSimilarity(S, 1, 100, 'same signal');
        S(:,:,2) = Y1(:,:,3);
        checkSimilarity(S, 1, 100, 'surrogate MVAR-gaussian');
        S(:,:,2) = Y2(:,:,3);
        checkSimilarity(S, 1, 100, 'surrogate MVAR-residuals');
        S(:,:,2) = Y3(:,:,3);
        checkSimilarity(S, 1, 100, 'surrogate mFT');
    end
end

function checkSimilarity(S, tRange, nDft, prefix)
    nodeNum = size(S,1);

    % calc mean difference (normalized)
    sm = squeeze(mean(S,2));
    dm = sm(:,1)-sm(:,2);
    TTp = abs(dm) / (tRange/2); % normalize

    % calc std deviation difference (normalized)
    ss = squeeze(std(S,1,2));
    ds = ss(:,1)-ss(:,2);
    FTp = abs(ds) / (tRange/4); % normalize

    % calc amplitude difference
    D = [];
    for j=1:2
        figure; D(:,:,j) = plotDft(S(:,:,j),nDft); % range [0, 1] signal
    end
    D = D * (nDft/2) / (tRange/2); % normalize
    E = D(:,:,1)-D(:,:,2);
    freqCos = getCosSimilarity(D(:,:,1),D(:,:,2));
    
    % calc auto-correlation
    D = [];
    for j=1:2
        figure; D(:,:,j) = plotAutoCorrelation(S(:,:,j),15);
    end
    acCos = getCosSimilarity(D(:,2:end,1),D(:,2:end,2));

    % calc FC
    nanx = eye(nodeNum,nodeNum);
    nanx(nanx==1) = NaN;
    for j=1:2
        FC(:,:,j) = calcFunctionalConnectivity(S(:,:,j));
    end
    cosSimFC = getCosSimilarity(FC(:,:,1) + nanx, FC(:,:,2) + nanx);

    % full correlation
    S1 = S(:,:,1); S2 = S(:,:,2);
    fullCorr = corr(S1(:),S2(:));

    % calc MSC
    for j=1:2
        MSC(:,:,:,j) = calcMSCoherence(S(:,:,j), [], [], [], 20);
    end
    cosSimMSC = getCosSimilarity(MSC(:,:,2:end-1,1) + nanx, MSC(:,:,2:end-1,2) + nanx);

    % calc PMSC
    for j=1:2
        PMSC(:,:,:,j) = calcPartialMSCoherence(S(:,:,j), [], [], [], 20);
    end
    cosSimPMSC = getCosSimilarity(PMSC(:,:,2:end-1,1) + nanx, PMSC(:,:,2:end-1,2) + nanx);

    % calc pGC
    for j=1:2
        GC(:,:,j) = calcPairwiseGCI(S(:,:,j), [], [], [], 4);
    end
    GC(isinf(GC) & GC>0) = 1.0e+50; % replace inf to dummy value
    GC(isinf(GC) & GC<0) = -1.0e+50; % replace inf to dummy value
    cosSimGC = getCosSimilarity(GC(:,:,1),GC(:,:,2));

    % calc mvGC
    for j=1:2
        mGC(:,:,j) = calcMultivariateGCI(S(:,:,j), [], [], [], 4);
    end
    cosSimMGC = getCosSimilarity(mGC(:,:,1),mGC(:,:,2));

    % calc NCC
    for j=1:2
        NCC(:,:,:,j) = calcCrossCorrelation(S(:,:,j), [], [], [], 4);
    end
    NCC(:,:,5,:) = [];
    cosSimNCC = getCosSimilarity(NCC(:,:,:,1) + nanx, NCC(:,:,:,2) + nanx);

    % calc NPCC
    for j=1:2
        NPCC(:,:,:,j) = calcPartialCrossCorrelation(S(:,:,j), [], [], [], 4);
    end
    NPCC(:,:,5,:) = [];
    cosSimNPCC = getCosSimilarity(NPCC(:,:,:,1) + nanx, NPCC(:,:,:,2) + nanx);

    disp([prefix ' : m=' num2str(mean(TTp)) ', s=' num2str(mean(FTp)) ', amp=' num2str(freqCos) ', ac=' num2str(acCos) ', fcorr=' num2str(fullCorr) ', cosFC=' num2str(cosSimFC) ', cosNCC=' num2str(cosSimNCC) ', cosNPCC=' num2str(cosSimNPCC) ', cosMSC=' num2str(cosSimMSC) ', cosPMSC=' num2str(cosSimPMSC) ', cosPGC=' num2str(cosSimGC) ', cosMGC=' num2str(cosSimMGC)]);
end


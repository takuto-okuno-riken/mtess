
function testSurrogateMVARCell
    % load signals
    load('data/ww32-1.mat');

    %% test pattern 1 
    for i=1:8, CX{i}=X; end
    for lags=1:5
        netMVAR = initMvarNetworkWithCell(CX, [], [], [], lags, 0.01);
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
    
    nanx = eye(nodeNum,nodeNum);
    nanx(nanx==1) = NaN;
    for j=1:2
        FC(:,:,j) = calcFunctionalConnectivity(S(:,:,j));
    end
    cosSimFC = getCosSimilarity(FC(:,:,1) + nanx, FC(:,:,2) + nanx);

    for j=1:2
        GC(:,:,j) = calcPairwiseGCI(S(:,:,j), [], [], [], 4);
    end
    GC(isinf(GC) & GC>0) = 1.0e+50; % replace inf to dummy value
    GC(isinf(GC) & GC<0) = -1.0e+50; % replace inf to dummy value
    cosSimGC = getCosSimilarity(GC(:,:,1),GC(:,:,2));

    disp([prefix ' : m=' num2str(mean(TTp)) ', s=' num2str(mean(FTp)) ', amp=' num2str(freqCos) ', cosFC=' num2str(cosSimFC) ', cosPGC=' num2str(cosSimGC)]);
end


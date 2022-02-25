% Before using this function, download ARR surrogates codes from
% https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/fMRI_dynamics/Liegeois2017_Surrogates
% and add a path "Liegeois2017_Surrogates" and sub folders.

function testSurrogateRL2017
    % load signals
    load('data/ww32-1.mat');

    %% test pattern 1 
    Y1 = CBIG_RL2017_get_AR_surrogate(X.',1,1,'gaussian');
    Y2 = CBIG_RL2017_get_AR_surrogate(X.',1,1,'nongaussian');
    Y3 = CBIG_RL2017_get_PR_surrogate(X.',1);
    Y1 = Y1.';
    Y2 = Y2.';
    Y3 = Y3.';
    Y4 = surrogateMFT(X);
    netMVAR = initMvarNetwork(X, [], [], [], 1);
    Y5 = surrogateMVAR(X, [], [], [], netMVAR);

    X1 = X(:,1:end-1);
    figure; plotTwoSignals(X,Y1);
    figure; plotTwoSignals(X,Y2);
    figure; plotTwoSignals(X1,Y3);
    figure; plotTwoSignals(X,Y4);
    figure; plotTwoSignals(X,Y5);

    % similarity check
    S = X;
    S(:,:,2) = X;
    checkSimilarity(S, 1, 100, 'same signal');
    S(:,:,2) = Y1;
    checkSimilarity(S, 1, 100, 'surrogate LS-ARR-gauss');
    S(:,:,2) = Y2;
    checkSimilarity(S, 1, 100, 'surrogate LS-ARR-nongauss');
    X1(:,:,2) = Y3;
    checkSimilarity(X1, 1, 100, 'surrogate LS-PR');
    S(:,:,2) = Y4;
    checkSimilarity(S, 1, 100, 'surrogate mFT');
    S(:,:,2) = Y5;
    checkSimilarity(S, 1, 100, 'surrogate MVAR');
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


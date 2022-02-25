
function testSurrogateAAFT
    % load signals
    load('data/ww32-1.mat');

    %% test pattern 1 
    Y1 = surrogateFT(X,10);
    Y2 = surrogateAAFT(X,10);
    Y3 = surrogateMFT(X,10);
    Y4 = surrogateMAAFT(X,10);
    Y5 = surrogateIAAFT(X,100,10);
    Y6 = surrogateMIAAFT(X,100,10);
    Y7 = surrogateRS(X,10);
    Y8 = surrogateMRS(X,10);
    Y9 = surrogateRG(X,10);
    Y10 = surrogateMRG(X,10);
    
    figure; plotTwoSignals(X,Y1(:,:,3));
    figure; plotTwoSignals(X,Y2(:,:,3));
    figure; plotTwoSignals(X,Y3(:,:,3));
    figure; plotTwoSignals(X,Y4(:,:,3));
    figure; plotTwoSignals(X,Y5(:,:,3));
    figure; plotTwoSignals(X,Y6(:,:,3));
    figure; plotTwoSignals(X,Y7(:,:,3));
    figure; plotTwoSignals(X,Y8(:,:,3));
    figure; plotTwoSignals(X,Y9(:,:,3));
    figure; plotTwoSignals(X,Y10(:,:,3));

    if Y1(:,:,3)==Y1(:,:,2), disp('bad surrogate : Y1'); end
    if Y2(:,:,3)==Y2(:,:,2), disp('bad surrogate : Y2'); end
    if Y3(:,:,3)==Y3(:,:,2), disp('bad surrogate : Y3'); end
    if Y4(:,:,3)==Y4(:,:,2), disp('bad surrogate : Y4'); end
    if Y5(:,:,3)==Y5(:,:,2), disp('bad surrogate : Y5'); end
    if Y6(:,:,3)==Y6(:,:,2), disp('bad surrogate : Y6'); end
    if Y7(:,:,3)==Y7(:,:,2), disp('bad surrogate : Y7'); end
    if Y8(:,:,3)==Y8(:,:,2), disp('bad surrogate : Y8'); end
    if Y9(:,:,3)==Y9(:,:,2), disp('bad surrogate : Y9'); end
    if Y10(:,:,3)==Y10(:,:,2), disp('bad surrogate : Y10'); end
    
    % similarity check
    S = X;
    S(:,:,2) = X;
    checkSimilarity(S, 1, 100, 'same signal');
    S(:,:,2) = Y9(:,:,3);
    checkSimilarity(S, 1, 100, 'surrogate RG');
    S(:,:,2) = Y7(:,:,3);
    checkSimilarity(S, 1, 100, 'surrogate RS');
    S(:,:,2) = Y1(:,:,3);
    checkSimilarity(S, 1, 100, 'surrogate FT');
    S(:,:,2) = Y2(:,:,3);
    checkSimilarity(S, 1, 100, 'surrogate AAFT');
    S(:,:,2) = Y5(:,:,3);
    checkSimilarity(S, 1, 100, 'surrogate IAAFT');
    S(:,:,2) = Y10(:,:,3);
    checkSimilarity(S, 1, 100, 'surrogate mRG');
    S(:,:,2) = Y8(:,:,3);
    checkSimilarity(S, 1, 100, 'surrogate mRS');
    S(:,:,2) = Y3(:,:,3);
    checkSimilarity(S, 1, 100, 'surrogate mFT');
    S(:,:,2) = Y4(:,:,3);
    checkSimilarity(S, 1, 100, 'surrogate mAAFT');
    S(:,:,2) = Y6(:,:,3);
    checkSimilarity(S, 1, 100, 'surrogate mIAAFT');
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
%{
    % show mean amplitudes
    E = nanmean([D(:,:,1);D(:,:,2)],1); f = 0.5*(0:(nDft/2))/nDft;
    figure; plot(f(2:end-1), E.');
    title('mean Single-Sided Amplitude Spectrum of S(t)'); xlabel('f (Hz)'); ylabel('mean |P1(f)|');
%}
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
    freqRange = 2:10;
    for j=1:2
        MSC(:,:,:,j) = calcMSCoherence(S(:,:,j), [], [], [], nDft);
    end
    cosSimMSC = getCosSimilarity(MSC(:,:,2:end-1,1) + nanx, MSC(:,:,2:end-1,2) + nanx);

    % calc PMSC
    for j=1:2
        PMSC(:,:,:,j) = calcPartialMSCoherence(S(:,:,j), [], [], [], nDft);
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


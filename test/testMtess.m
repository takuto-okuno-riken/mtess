function testMtess
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    uuOrg = si;
    nodeNum = 8;
    sigLen = 500;

    % test with net-pat3-8x500-idx6-* data (DCM simulated fMRI signal)
    load('data/cx-8x500-idx6.mat');

    range = [0 1];
    acLags = 5;
    pacLags = 13;
    ccLags = 2;
    pccLags = 4;
    pccFunc = @calcPartialCrossCorrelation;
    [MTS, MTSp, nMTS, nMTSp, Means, Stds, ACs, PACs, FCs, PCs, CCs, PCCs] = calcMtess(CX, range, pccFunc, acLags, pacLags, ccLags, pccLags);
    plotMtessAllMatrix(MTS, MTSp, 'cx-8x500-idx6');
    P=squeeze(MTSp(1,2:8,:));
    figure; plotMtessSpiderPlot(P);
    figure; h=bar3(MTS); 
    % cache version
    CXNames = {};
    for i=1:length(CX), CXNames{i} = ['dmy' num2str(i)]; end
    [MTS2, MTSp2, nMTS2, nMTSp2] = calcMtess_c(CX, range, pccFunc, acLags, pacLags, ccLags, pccLags, CXNames);
    [MTS3, MTSp3, nMTS3, nMTSp3] = calcMtessCross_c(CX, CX, range, pccFunc, acLags, pacLags, ccLags, pccLags, CXNames, CXNames);

%    figure; plotMtessPolarplot(MTSp, 1);
    P=squeeze(MTSp(1,2:8,:));
    figure; plotMtessSpiderPlot(P);
    title('cx-8x500-idx6 : MTESS polar chart');

    % 3D bar graph
    figure; h=bar3(MTS); title('MTESS matrix 3D bar graph');
    for i=1:length(h) h(i).FaceAlpha=0.6; h(i).EdgeAlpha=0.6; end
    xlabel('Cell number');
    ylabel('Cell number');
    zlabel('MTESS'); zticks([0 1 2 3 4 5]);

    % check 'half' data input
    for i=1:length(CX), CX{i}=half(CX{i}); end
    [MTS, MTSp, nMTS, nMTSp, Means, Stds, ACs, PACs, FCs, PCs, CCs, PCCs] = calcMtess(CX, range, pccFunc, acLags, pacLags, ccLags, pccLags);
    P=squeeze(MTSp(1,2:8,:));
    figure; plotMtessSpiderPlot(P);
    title('cx-8x500-idx6(half) : MTESS polar chart');
        
    % test surrogate algorithms
    for i=1:length(CX), CX{i}=single(CX{i}); end
    algoNum = 7;
    names = {'original', 'mRG', 'mRS', 'mFT', 'mAAFT', 'VAR4', 'FT'};
    CX2{1} = CX{1};
    CX2{2} = single(surrogateMRG(CX{1}, 1));
    CX2{3} = single(surrogateMRS(CX{1}, 1));
    CX2{4} = single(surrogateMFT(CX{1}, 1));
    CX2{5} = single(surrogateMAAFT(CX{1}, 1));
    net = initMvarNetwork(CX{1}, [], [], [], 4);
    CX2{6} = single(surrogateMVAR(CX{1}, [], [], [], net, 'gaussian', 1));
    CX2{7} = single(surrogateFT(CX{1}, 1));
%    CXt = CX; CX = CX2;
%    save('data/cx-8x500-idx6-surrogate.mat','CX','names');
    [MTS, MTSp, nMTS, nMTSp, Means, Stds, ACs, PACs, FCs, PCs, CCs, PCCs] = calcMtess(CX2, range, pccFunc, acLags, pacLags, ccLags, pccLags);
    plotMtessAllMatrix(MTS, MTSp, 'surrogate algorithms');
    A=squeeze(nMTS(1,2:algoNum,:));
    figure; plot(A.'); title('surrogate algorithms : node MTESS');
    legend(names(2:algoNum)); ylim([0,5]);

%    figure; plotMtessPolarplot(MTSp, 1);
    P=squeeze(MTSp(1,2:algoNum,:));
    figure; plotMtessSpiderPlot(P);
    legend(names(2:algoNum), 'Location', 'northeast');
    title('surrogate algorithms : MTESS polar chart');

    for i=2:algoNum
        P=squeeze(nMTSp(1,i,:,:));
        figure; plotMtessSpiderPlot(P);
        title(['surrogate algorithms ' names{i} ' : node MTESS polar chart']);
    end

    % synthetic lines
    CX3{1} = CX{1};
    CX3{2} = single(ones(nodeNum,sigLen)) * 0.5; % flat @ 0.5
    CX3{3} = single(ones(nodeNum,sigLen)) * 0.9; % flat @ 0.9
    CX3{4} = single(uuOrg(1:nodeNum,1:sigLen)); % random
    y = sin([1:sigLen] * pi / 8) * 0.5 + 0.5;
    CX3{5} = single(repmat(y, [nodeNum, 1]));
    figure; plotTwoSignals(CX3{1},CX3{2},false,[0,1]);
%{
    % check flat line
    figure; FC = calcFunctionalConnectivity(CX3{2});
    figure; PC2 = calcPartialCorrelation__(CX3{2});
    % check sin curve
    figure; FC = plotFunctionalConnectivity(CX3{5});
    figure; PC = plotPartialCorrelation(CX3{5});
    figure; PC2 = calcPartialCorrelation__(CX3{5});
    PCC = calcPartialCrossCorrelation(CX3{5}, [], [], [], pccLags);
    PC3 = PCC(:,:,pccLags+1);
%}
    [MTS, MTSp, nMTS, nMTSp, Means, Stds, ACs, PACs, FCs, PCs, CCs, PCCs] = calcMtess(CX3, range, pccFunc, acLags, pacLags, ccLags, pccLags);
    P=squeeze(MTSp(1,2:5,:));
    figure; plotMtessSpiderPlot(P);
    legend('flat@0.5', 'flat@0.9', 'random', 'sin');
    title('synthetic lines : MTESS polar chart');

    % node effect checking
    CX3 = CX2;
    CX3{4}(4,:)=0.5; % flat @ 0.5
    CX3{4}(5,:)=0.9; % flat @ 0.9
    CX3{4}(6,:)=uuOrg(6,1:sigLen); % random
    CX3{4}(7,:)=sin([1:sigLen] * pi / 8) * 0.5 + 0.5; % sin
    figure; plotTwoSignals(CX2{4},CX3{4},false,[0,1]);
    CX3{5}(4,:)=0.5; % flat @ 0.5
    CX3{5}(5,:)=0.9; % flat @ 0.9
    CX3{5}(6,:)=uuOrg(6,1:sigLen); % random
    CX3{5}(7,:)=sin([1:sigLen] * pi / 8) * 0.5 + 0.5; % sin
    figure; plotTwoSignals(CX2{4},CX3{4},false,[0,1]);

    [MTS, MTSp, nMTS, nMTSp, Means, Stds, ACs, PACs, FCs, PCs, CCs, PCCs] = calcMtess(CX3, range, pccFunc, acLags, pacLags, ccLags, pccLags);
    plotMtessAllMatrix(MTS, MTSp, 'node effect checking');
    A=squeeze(nMTS(1,2:6,:));
    figure; plot(A.'); title('node effect checking : node MTESS');
    legend('mRG', 'mRS', 'mFT_mod', 'mAAFT_mod', 'VAR4'); ylim([0,5]);

    P=squeeze(MTSp(1,2:6,:));
    figure; plotMtessSpiderPlot(P);
    legend('mRG', 'mRS', 'mFT_mod', 'mAAFT_mod', 'VAR4', 'Location', 'northeast');
    title('node effect checking : MTESS polar chart');
    
    for i=4:5
        P=squeeze(nMTSp(1,i,:,:));
        figure; plotMtessSpiderPlot(P);
        title(['node effect checking ' names{i} ' : node MTESS polar chart']); legend;
    end

    % test with ww32-* data (TVB simulated neural signal)
    CX = {};
    for i=1:4
        load(['data/ww32-' num2str(i) '.mat']);
        figure; plot(X.'); title(['ww32-' num2str(i)]);
        CX{i} = X;
    end

    [MTS, MTSp, nMTS, nMTSp, Means, Stds, ACs, PACs, FCs, PCs, CCs, PCCs] = calcMtess(CX, range, pccFunc, acLags, pacLags, ccLags, pccLags);
    plotMtessAllMatrix(MTS, MTSp, 'ww32-* data');

%    figure; plotMtessPolarplot(MTSp, 1);
    P=squeeze(MTSp(1,2:4,:));
    figure; plotMtessSpiderPlot(P);
    
    % real HCP fMRI signal
%{
    signals = connData2signalsFile('', {}, 'me', 'data/hcp', 'hcp');
    CX = {};
    for i=1:6
        CX{i} = convert2SigmoidSignal(single(signals{i}));
    end
    save('data/hcp-metest-132x1190.mat','CX');
%}
    filename = 'results/hcp-metest-132x1190-result.mat';
    if exist(filename,'file')
        load(filename);
    else
        load('data/hcp-metest-132x1190.mat');
        pccLags = 2;
        pccFunc = @calcSvPartialCrossCorrelation;
        [MTS, MTSp, nMTS, nMTSp, Means, Stds, ACs, PACs, FCs, PCs, CCs, PCCs] = calcMtess(CX, range, pccFunc, acLags, pacLags, ccLags, pccLags);
        save(filename,'MTS', 'MTSp', 'nMTS', 'nMTSp', 'Means', 'Stds', 'ACs', 'PACs', 'FCs', 'PCs', 'CCs', 'PCCs');
    end
    plotMtessAllMatrix(MTS, MTSp, 'real HCP fMRI signal');
    P=squeeze(MTSp(1,2:6,:));
    figure; plotMtessSpiderPlot(P); title('real HCP fMRI : MTESS polar chart');
    
    A=squeeze(nMTS(1,2:6,:));
    figure; plot(A.'); title('real HCP fMRI : node MTESS');
end

function plotMtessPolarplot(MTSp, src)
    % male vs. female MTESS polarplot
    theta = 0:(2*pi)/7:2*pi;
    C = [];
    for j=2:size(MTSp,2)
        for i=1:7, A=MTSp(src,j,i); C(i,j)=nanmean(A(:)); end
    end
    C(8,:)=C(1,:);
    polarplot(theta,C(:,1));
    for j=2:size(MTSp,2)
        hold on; polarplot(theta,C(:,j)); hold off;
    end
    title(['MTESS radar chart']); legend;
end

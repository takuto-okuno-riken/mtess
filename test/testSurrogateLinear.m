%
% None-Linear test by J.Theilear and D.Prichard, Physica D: Nonlinear Phenomena (1996) pp.221-235.
% using FT surrogate and etc.
%

function testSurrogateLinear
    % load signals
    load('data/ww32-1.mat');

    %% univariate Linear test by FT surrogate
    Y = surrogateFT(X,399);
%    figure; plotTwoSignals(X,Y(:,:,1));

    % test linearity
    statisticParams = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    for i=1:32
        figure; [H(i), P(i), T(i,:), Rank(i)] = plotSurrogateTest(X(i,:), Y(i,:,:), @calcSurrStatLinear, statisticParams);
    end
    [H2, P2, T2, Rank2] = calcSurrogateTest(X, Y, @calcSurrStatLinear, statisticParams);
    if Rank(:) == Rank2
        disp('got same result each vs. matrix');
    else
        disp('did not get same result each vs. matrix');
    end
    disp(['significantly not linear (' num2str(sum(H2)) ' / 32) by FT surrogate']);


    %% univariate Linear test by AAFT surrogate
    Y = surrogateAAFT(X,399);
    [H2, P2, T2, Rank2] = calcSurrogateTest(X, Y, @calcSurrStatLinear, statisticParams);
    disp(['significantly not linear (' num2str(sum(H2)) ' / 32) by AAFT surrogate']);


    %% univariate Linear test by AR(6) surrogate
    Y = surrogateAR(X,6,'gaussian',399);
    figure; plotTwoSignals(X,Y(:,:,1));
    for i=1:32
        figure; [H(i), P(i), T(i,:), Rank(i)] = plotSurrogateTest(X(i,:), Y(i,:,:), @calcSurrStatLinear, statisticParams);
    end
    disp(['significantly not linear (' num2str(sum(H)) ' / 32) by AR(6) surrogate']);


    %% univariate Linear test by mFT surrogate
    Y = surrogateMFT(X,399);
    [H2, P2, T2, Rank2] = calcSurrogateTest(X, Y, @calcSurrStatLinear, statisticParams);
    disp(['significantly not linear (' num2str(sum(H2)) ' / 32) by mFT surrogate']);

    %% univariate Linear test by mAAFT surrogate
    Y = surrogateMAAFT(X,399);
    [H2, P2, T2, Rank2] = calcSurrogateTest(X, Y, @calcSurrStatLinear, statisticParams);
    disp(['significantly not linear (' num2str(sum(H2)) ' / 32) by mAAFT surrogate']);

    %% univariate Linear test by VAR(3) surrogate
    netMVAR = initMvarNetwork(X, [], [], [], 3);
    Y = surrogateMVAR(X, [], [], [], netMVAR, 'gaussian', 399);
    [H2, P2, T2, Rank2] = calcSurrogateTest(X, Y, @calcSurrStatLinear, statisticParams);
    disp(['significantly not linear (' num2str(sum(H2)) ' / 32) by VAR(3) surrogate']);
end


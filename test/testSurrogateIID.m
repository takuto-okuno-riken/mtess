%
% IID (independent and identically distributed) test by J.Theilear and D.Prichard, Physica D: Nonlinear Phenomena (1996) pp.221-235.
% using random shuffling surrogate and etc.
%

function testSurrogateIID
    % load signals
    load('data/ww32-1.mat');

    %% univariate IID test by RS surrogate
    Y = surrogateRS(X,399);
    figure; plotTwoSignals(X,Y(:,:,1));

    % test none linearity
    for i=1:32
        figure; [H(i), P(i), T(i,:), Rank(i)] = plotSurrogateTest(X(i,:), Y(i,:,:), @calcSurrStatIID, []);
    end
    [H2, P2, T2, Rank2] = calcSurrogateTest(X, Y, @calcSurrStatIID, []);
    if Rank(:) == Rank2
        disp('got same result each vs. matrix');
    else
        disp('did not get same result each vs. matrix');
    end
    disp(['significantly not IID (' num2str(sum(H2)) ' / 32) by RS surrogate']);


    %% check just random gaussian
    X = normrnd(0, 1, 32, 300);
    Y = surrogateRS(X,399);
    [H2, P2, T2, Rank2] = calcSurrogateTest(X, Y, @calcSurrStatIID, []);
    disp(['significantly not IID (' num2str(sum(H2)) ' / 32) by RS surrogate']);
end

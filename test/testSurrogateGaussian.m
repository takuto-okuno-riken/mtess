%
% Gaussian distribution test by J.Theilear and D.Prichard, Physica D: Nonlinear Phenomena (1996) pp.221-235.
% using random gaussian surrogate and etc.
%

function testSurrogateGaussian
    % load signals
    load('data/ww32-1.mat');

    %% univariate Gaussian test by RG surrogate
    Y = surrogateRG(X,399);
    figure; plotTwoSignals(X,Y(:,:,1));

    % test gaussian
    for i=1:32
        figure; [H(i), P(i), T(i,:), Rank(i)] = plotSurrogateTest(X(i,:), Y(i,:,:), @calcSurrStatGaussian, []);
    end
    [H2, P2, T2, Rank2] = calcSurrogateTest(X, Y, @calcSurrStatGaussian, []);
    if Rank(:) == Rank2
        disp('got same result each vs. matrix');
    else
        disp('did not get same result each vs. matrix');
    end
    disp(['significantly not gaussian (' num2str(sum(H2)) ' / 32) by RG surrogate']);

    %% check just random gaussian
    X = normrnd(0, 1, 32, 300);
    Y = surrogateRG(X,399);
    [H2, P2, T2, Rank2] = calcSurrogateTest(X, Y, @calcSurrStatIID, []);
    disp(['significantly not gaussian (' num2str(sum(H2)) ' / 32) by RS surrogate']);
end

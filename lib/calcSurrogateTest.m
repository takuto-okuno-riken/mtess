%%
% Caluclate surrogate rank test
% returns significance vector (H=1 or 0), p-value vector (P), discriminating statistic matrix T (node x surrNum) and original rank vector(Rank).
% input:
%  X                original multivariate time series matrix (node x time series)
%  Y                surrogate multivariate time series matrix (node x time series x surrNum)
%  statisticFunc    discriminating statistic function
%  statisticParams  discriminating statistic function parameters
%  side             bottm-side(1), both-side(2), top-side(3) (default:2)
%  alpha            the significance level of statistic (default:0.05)

function [H, P, T, Rank] = calcSurrogateTest(X, Y, statisticFunc, statisticParams, side, alpha)
    if nargin < 6, alpha = 0.05; end
    if nargin < 5, side = 2; end
    if nargin < 4, statisticParams = []; end
    nodeNum = size(Y,1);
    surrNum = size(Y,3);
    T = nan(nodeNum,surrNum);
    
    % calculate discriminating statistic
    for i=1:nodeNum
        T(i,1) = statisticFunc(X(i,:), statisticParams);
        for k=1:surrNum
            T(i,k+1) = statisticFunc(Y(i,:,k), statisticParams);
        end
    end
    
    % rank test
    [B, I] = sort(T,2);
    R = repmat([1:surrNum+1],nodeNum,1);
    I(I~=1) = 0;
    R2 = R .* I;
    Rank = sum(R2,2);
    if side == 1
        P = Rank / (surrNum + 1);
    else if side == 3
        P = 1 - ((Rank-1) / (surrNum + 1));
    else
        R3 = Rank;
        n = (surrNum + 1) / 2;
        idx = find(R3 > n);
        R3(idx) = (surrNum + 2) - R3(idx);
        P = R3 / n;
    end
    H = P < alpha;
end

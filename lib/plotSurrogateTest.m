%%
% Plot surrogate rank test
% returns significance vector (H=1 or 0), p-value vector (P), discriminating statistic matrix T (node x surrNum) and original rank vector(Rank).
% input:
%  X                original multivariate time series matrix (node x time series)
%  Y                surrogate multivariate time series matrix (node x time series x surrNum)
%  statisticFunc    discriminating statistic function
%  statisticParams  discriminating statistic function parameters
%  side             bottm-side(1), both-side(2), top-side(1) (default:2)
%  alpha            the significance level of statistic (default:0.05)

function [H, P, T, Rank] = plotSurrogateTest(X, Y, statisticFunc, statisticParams, side, alpha)
    if nargin < 6, alpha = 0.05; end
    if nargin < 5, side = 2; end
    if nargin < 4, statisticParams = []; end
    
    [H, P, T, Rank] = calcSurrogateTest(X, Y, statisticFunc, statisticParams, side, alpha);
    
    nodeNum = size(Y,1);
    surrNum = size(Y,3);
    for i=1:nodeNum
        for k=2:surrNum
            hold on;
            plot([T(i,k),T(i,k)],[i-1,i-0.5],'blue');
            hold off;
        end
        hold on;
        plot([T(i,1),T(i,1)],[i-1,i],'red');
        hold off;
    end
end

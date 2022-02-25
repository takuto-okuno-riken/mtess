%%
% Surrogate univariate signal generation by Random Shuffling (RS)
% input and output is multivariate, but it is processed by univariate.
% returns surrogated signals Y (node x time seriese x surrNum)
% input:
%  X          multivariate time series matrix (node x time series)
%  surrNum    output number of surrogate samples (default:1)

function Y = surrogateRS(X, surrNum)
    if nargin < 2, surrNum = 1; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        for i=1:nodeNum
            x = X(i,:);
            perm = randperm(sigLen);
            Y(i,:,k) = x(perm);
        end
    end
end

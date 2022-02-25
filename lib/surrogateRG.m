%%
% Surrogate univariate signal generation by Random Gaussian (RG)
% input and output is multivariate, but it is processed by univariate.
% returns surrogated signals Y (node x time seriese x surrNum)
% input:
%  X          multivariate time series matrix (node x time series)
%  surrNum    output number of surrogate samples (default:1)

function Y = surrogateRG(X, surrNum)
    if nargin < 2, surrNum = 1; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    Y = nan(nodeNum,sigLen,surrNum);
    for i=1:nodeNum
        m = mean(X(i,:));
        s = std(X(i,:),1);
        for k=1:surrNum
            Y(i,:,k) = normrnd(m,s,1,sigLen);
        end
    end
end

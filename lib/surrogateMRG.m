%%
% Surrogate univariate signal generation by Random Gaussian (RG)
% returns surrogated signals Y (node x time seriese x surrNum)
% input:
%  X          multivariate time series matrix (node x time series)
%  surrNum    output number of surrogate samples (default:1)

function Y = surrogateMRG(X, surrNum)
    if nargin < 2, surrNum = 1; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    M  = mean(X.');
    EC = cov(X.',1);
    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        Y(:,:,k) = (mvnrnd(M,EC,size(X,2))).';
    end
end

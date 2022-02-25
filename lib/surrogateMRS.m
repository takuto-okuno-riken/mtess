%%
% Surrogate multivariate signal generation by Random Shuffling (RS)
% returns surrogated signals Y (node x time seriese x surrNum)
% input:
%  X          multivariate time series matrix (node x time series)
%  surrNum    output number of surrogate samples (default:1)

function Y = surrogateMRS(X, surrNum)
    if nargin < 2, surrNum = 1; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        perm = randperm(sigLen);
        for i=1:nodeNum
            x = X(i,:);
            Y(i,:,k) = x(perm);
        end
    end
end

%%
% Surrogate univariate signal generation by Amplitude adjusted Fourier Transform (AAFT)
% input and output is multivariate, but it is processed by univariate.
% returns surrogated signals (Y)
% input:
%  X          multivariate time series matrix (node x time series)
%  surrNum    output number of surrogate samples (default:1)

function Y = surrogateAAFT(X, surrNum)
    if nargin < 2, surrNum = 1; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        for i=1:nodeNum
            Y(i,:,k) = surrogateAAFTuni(X(i,:));
        end
    end
end

function y = surrogateAAFTuni(x)
    n = length(x);
    [x2, I] = sort(x);
    [I2, J] = sort(I);
    r = sort(randn(1,n));
    g = r(J);
    h = surrogateFT(g);
    [h2, K] = sort(h);
    [K2, L] = sort(K);
    y = x2(L);
end

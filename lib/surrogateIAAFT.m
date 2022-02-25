%%
% Surrogate univariate signal generation by Iterated Amplitude adjusted Fourier Transform (IAAFT)
% input and output is multivariate, but it is processed by univariate.
% returns surrogated signals (Y)
% input:
%  X          multivariate time series matrix (node x time series)
%  maxIter    maximum iteration number (default:100)
%  surrNum    output number of surrogate samples (default:1)

function Y = surrogateIAAFT(X, maxIter, surrNum)
    if nargin < 3, surrNum = 1; end
    if nargin < 2, maxIter = 100; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        for i=1:nodeNum
            Y(i,:,k) = surrogateIAAFTuni(X(i,:), maxIter);
        end
    end
end

function y = surrogateIAAFTuni(x, maxIter)
    [x2, I] = sort(x);
    S = abs(fft(x));
    y = surrogateAAFT(x);
    for i=1:maxIter
        R = fft(y);
        sm = (R ./ abs(R)) .* S;
        h = ifft(sm);
        [h2, K] = sort(h);
        [K2, L] = sort(K);
        y2 = x2(L);
        if y == y2, break; end
        y = y2;
    end
end

%%
% Surrogate univariate signal generation by Fourier Transform (FT)
% input and output is multivariate, but it is processed by univariate.
% returns surrogated signals Y (node x time seriese x surrNum)
% input:
%  X          multivariate time series matrix (node x time series)
%  surrNum    output number of surrogate samples (default:1)

function Y = surrogateFT(X, surrNum)
    if nargin < 2, surrNum = 1; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        for i=1:nodeNum
            Y(i,:,k) = surrogateFTuni(X(i,:));
        end
    end
end

function y = surrogateFTuni(x)
    n = length(x);
    if mod(n, 2) == 0
        r = exp(2i*pi*rand(1, n/2 - 1));
        v = [1, r, 1, flip(conj(r))];
    else
        r = exp(2i*pi*rand(1, (n-1)/2));
        v = [1, r, flip(conj(r))];
    end
    xf = fft(x);
    y = ifft(xf .* v);
end

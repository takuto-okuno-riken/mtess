%%
% Surrogate multivariate signal generation by Fourier Transform (FT)
% returns surrogated signals Y (node x time seriese x surrNum)
% D. Prichard and J. Theiler, Generating surrogate data for time series with several simultaneously measured variables, Phys. Rev. Lett. 73, 951 (1994). 
% input:
%  X          multivariate time series matrix (node x time series)
%  surrNum    output number of surrogate samples (default:1)

function Y = surrogateMFT(X, surrNum)
    if nargin < 2, surrNum = 1; end
    nodeNum = size(X,1);
    sigLen = size(X,2);

    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        if mod(sigLen, 2) == 0
            r = exp(2i*pi*rand(1, sigLen/2 - 1));
            v = [1, r, 1, flip(conj(r))];
        else
            r = exp(2i*pi*rand(1, (sigLen-1)/2));
            v = [1, r, flip(conj(r))];
        end
        Yk = fft(X,[],2);
        Y(:,:,k) = ifft(Yk .* repmat(v,nodeNum,1),[],2);
    end
end

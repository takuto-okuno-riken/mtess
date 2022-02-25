%%
% Surrogate multivariate signal generation by Iterated Amplitude adjusted Fourier Transform (AAFT)
% returns surrogated signals Y (node x time seriese x surrNum)
% input:
%  X          multivariate time series matrix (node x time series)
%  maxIter    maximum iteration number (default:100)
%  surrNum    output number of surrogate samples (default:1)

function Y = surrogateMIAAFT(X, maxIter, surrNum)
    if nargin < 3, surrNum = 1; end
    if nargin < 2, maxIter = 100; end

    nodeNum = size(X,1);
    [X2, I] = sort(X,2);
    S = abs(fft(X,[],2));
    Y = surrogateMAAFT(X, surrNum);
    for k=1:surrNum
        Yk = Y(:,:,k);
        for i=1:maxIter
            R = fft(Yk,[],2);
            sm = (R ./ abs(R)) .* S;
            H = ifft(sm,[],2);
            [H2, K] = sort(H,2);
            [K2, L] = sort(K,2);
            Y2 = [];
            for j=1:nodeNum
                t = X2(j,:);
                Y2(j,:) = t(L(j,:));
            end
            if Yk == Y2, break; end
            Yk = Y2;
        end
        Y(:,:,k) = Yk;
    end
end


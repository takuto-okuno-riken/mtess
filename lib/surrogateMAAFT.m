%%
% Surrogate multivariate signal generation by Amplitude adjusted Fourier Transform (AAFT)
% returns surrogated signals Y (node x time seriese x surrNum)
% input:
%  X          multivariate time series matrix (node x time series)
%  surrNum    output number of surrogate samples (default:1)

function Y = surrogateMAAFT(X, surrNum)
    if nargin < 2, surrNum = 1; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    
    [X2, I] = sort(X,2);
    [I2, J] = sort(I,2);
    R = sort(randn(nodeNum,sigLen),2);
    G = [];
    for i=1:nodeNum
        t = R(i,:);
        G(i,:) = t(J(i,:));
    end
    H = surrogateMFT(G, surrNum);
    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        [h2, K] = sort(H(:,:,k),2);
        [K2, L] = sort(K(:,:),2);
        for i=1:nodeNum
            t = X2(i,:);
            Y(i,:,k) = t(L(i,:));
        end
    end
end

%%
% Surrogate eigenvalue time-series (Shinn et al., 2022)
% returns surrogated signals Y (node x time seriese x surrNum)
% input:
%  X          multivariate time series matrix (node x time series)
%  surrNum    output number of surrogate samples (default:1)

function Y = surrogateEigenTS(X, surrNum)
    if nargin < 2, surrNum = 1; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    C = calcFunctionalConnectivity(X);
    S = surrogateEigenMatrix(C);
    Q = sqrtm(S);
    Y = nan(nodeNum,sigLen,surrNum);
    for i=1:surrNum
        Y(:,:,i) = Q * randn(nodeNum,sigLen);
    end
end

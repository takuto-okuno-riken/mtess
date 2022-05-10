%%
% Surrogate univariate signal generation by Auto-Regression (AR)
% input and output is multivariate, but it is processed by univariate.
% returns surrogated signals Y (node x time seriese x surrNum),
% AR coefficients (b) and AR residuals (r)
% input:
%  X          multivariate time series matrix (node x time series)
%  lags       number of lags for autoregression (default:3)
%  dist       distribution of noise to yield surrogates ('gaussian'(default), 'residuals')
%  surrNum    output number of surrogate samples (default:1)

function [Y, b, r] = surrogateAR(X, lags, dist, surrNum)
    if nargin < 4, surrNum = 1; end
    if nargin < 3, dist = 'gaussian'; end
    if nargin < 2, lags = 3; end

    nodeNum = size(X,1);
    sigLen = size(X,2);

    Y = nan(nodeNum,sigLen,surrNum);
    for i=1:nodeNum
        x = X(i,:);
        y = flipud(x.'); % need to flip signal
        for p=1:lags
            yj(:,p) = y(1+p:sigLen-lags+p);
        end
        % auto-regression (AR)
        xt = y(1:sigLen-lags);
        xti = [yj, ones(sigLen-lags,1)];
        % apply the regress function
        [b,r,~] = regressLinear(xt,xti);


        if strcmp(dist,'gaussian')
            m  = mean(r);
            sig = std(r,1);
            noise = normrnd(m,sig,length(r),1);
        else
            noise = r;
        end

        for k=1:surrNum
            S = x;
            perm      = randperm(sigLen-lags);
            S(1:lags) = S(perm(1):perm(1)+lags-1); %Initialization of the AR surrogate
            perm2     = randperm(sigLen-lags);

            for t=lags+1:sigLen
                S2 = [S(t-1:-1:t-lags), 1]; % might not be good to add bias

                % yield next time step
                S(t) = S2 * b + noise(perm2(t-lags));
            end
            Y(i,:,k) = S;
        end
    end
end

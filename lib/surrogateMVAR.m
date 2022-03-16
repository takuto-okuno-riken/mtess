%%
% Surrogate multivariate signal generation by multivariate VAR
% based on autoregressive (AR) surrogates (R. Liegeois et al., 2017)
% returns surrogated signals (node x time series x surrNum)(Y)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mVAR network
%  dist         distribution of noise to yield surrogates ('gaussian'(default), 'residuals')
%  surrNum      output number of surrogate samples (default:1)
%  yRange       range of Y value (default:[-0.2 1.2])

function Y = surrogateMVAR(X, exSignal, nodeControl, exControl, net, dist, surrNum, yRange)
    if nargin < 8, yRange = [-0.2 1.2]; end
    if nargin < 7, surrNum = 1; end
    if nargin < 6, dist = 'gaussian'; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    lags = net.lags;

    % set node input
    Xorg = [X; exSignal];

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    idxs = {};
    rvlen = length(net.rvec{1});
    for i=2:nodeNum, if rvlen > length(net.rvec{i}), rvlen = length(net.rvec{i}); end; end
    Err = single(nan(nodeNum,rvlen));
    for i=1:nodeNum
        for k=1:lags
            [~,idxs{i,k}] = find(control(i,:,k)==1);
        end
        Err(i,:) = net.rvec{i}(1:rvlen);
    end

    if strcmp(dist,'gaussian')
        P  = mean(Err.');
        EC = cov(Err.');
        noise = (mvnrnd(P,EC,size(Err,2)))';
    else
        noise = Err;
    end

    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        disp(['surrogate sample : ' num2str(k)]);
        S = single(Xorg);
        perm        = single(randperm(sigLen-lags));
        S(1:nodeNum,1:lags) = S(1:nodeNum,perm(1):perm(1)+lags-1); % Initialization of the AR surrogate
        perm2       = single(randperm(size(Err,2)));

        for t=lags+1:sigLen
            A = S(:,t);
            for i=1:nodeNum
    %        parfor i=1:nodeNum
                S2 = [];
                for p=1:lags
                    S2 = [S2; S(idxs{i,p},t-p)];
                end
                S2 = [S2; 1]; % might not be good to add bias

                % yield next time step
                A(i) = S2.' * net.bvec{i} + noise(i,perm2(t-lags));
            end
            % fixed over shoot values
            A(A < yRange(1)) = yRange(1);
            A(A > yRange(2)) = yRange(2);
            S(:,t) = A;
        end
        Y(:,:,k) = S(1:nodeNum,:);
    end
end

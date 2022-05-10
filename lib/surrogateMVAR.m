%%
% Surrogate multivariate signal generation by multivariate VAR
% based on autoregressive (AR) surrogates (R. Liegeois et al., 2017)
% returns surrogated signals (Y)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mVAR network
%  dist         distribution of noise to yield surrogates ('gaussian'(default), 'residuals')
%  surrNum      output number of surrogate samples (default:1)
%  yRange       range of Y value (default:[Xmin-Xrange/5, Xmax+Xrange/5])
%  Cin          coefficient matrix of VAR (default:[])
%  usegpu       use gpu for regress function (default:false)

function [Y, C] = surrogateMVAR(X, exSignal, nodeControl, exControl, net, dist, surrNum, yRange, Cin, usegpu)
    if nargin < 10, usegpu = false; end
    if nargin < 9, Cin = []; end
    if nargin < 8, yRange = NaN; end
    if nargin < 7, surrNum = 1; end
    if nargin < 6, dist = 'gaussian'; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    inputNum = nodeNum + exNum;
    lags = net.lags;

    % set node input
    Xorg = [X; exSignal];

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    % calc Y range
    if isnan(yRange)
        t = max(X(:)); d = min(X(:));
        r = t - d;
        yRange = [d-r/5, t+r/5];
    end

    rvlen = length(net.rvec{1});
    for i=2:nodeNum, if rvlen > length(net.rvec{i}), rvlen = length(net.rvec{i}); end; end
    Err = single(nan(nodeNum,rvlen));
    for i=1:nodeNum
        Err(i,:) = net.rvec{i}(1:rvlen);
    end
    % get coefficient matrix
    if isempty(Cin)
        C = single(zeros(nodeNum,inputNum*lags+1));
        for i=1:nodeNum
            idx = find(control(i,:,:)==1);
            C(i,[idx(:).' end]) = net.bvec{i};
        end
    else
        C = Cin;
    end

    if strcmp(dist,'gaussian')
        P  = mean(Err.');
        EC = cov(Err.');
        noise = (mvnrnd(P,EC,size(Err,2)))';
    else
        noise = Err;
    end
    S2 = ones(inputNum*lags+1,1);
    % use gpu array
    if usegpu
        C = gpuArray(C);
        S2 = gpuArray(S2);
        noise = gpuArray(noise);
    end

    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        disp(['surrogate sample : ' num2str(k)]);
        S = single(Xorg);
        perm        = single(randperm(sigLen-lags));
        S(1:nodeNum,1:lags) = S(1:nodeNum,perm(1):perm(1)+lags-1); % Initialization of the AR surrogate
        perm2       = single(randperm(size(Err,2)));

        for t=lags+1:sigLen
            A = S(:,t); % next output

            for p=1:lags
                S2(1+inputNum*(p-1):inputNum*p) = S(:,t-p);
            end
            A(1:nodeNum) = C * S2 + noise(:,perm2(t-lags));

            % fixed over shoot values
            if ~isempty(yRange)
                A(A < yRange(1)) = yRange(1);
                A(A > yRange(2)) = yRange(2);
            end
            S(:,t) = A;
        end
        Y(:,:,k) = S(1:nodeNum,:);
    end
end

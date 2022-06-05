%%
% Surrogate multivariate signal generation by multivariate PCVAR
% returns surrogated signals (Y)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mPCVAR network
%  dist         distribution of noise to yield surrogates ('gaussian'(default), 'residuals')
%  surrNum      output number of surrogate samples (default:1)
%  yRange       range of Y value (default:[Xmin-Xrange/5, Xmax+Xrange/5])

function Y = surrogateMpcvar(X, exSignal, nodeControl, exControl, net, dist, surrNum, yRange)
    if nargin < 8, yRange = NaN; end
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

    % calc Y range
    if isnan(yRange)
        t = max(X(:)); d = min(X(:));
        r = t - d;
        yRange = [d-r/5, t+r/5];
    end

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
        EC = cov(Err.',1);
        noise = (mvnrnd(P,EC,size(Err,2)))';
    else
        noise = Err;
    end

    mc = net.maxComp;
    mu = net.mu;
    if size(net.coeff{1},1) == size(net.coeff{1},2)
        for i=1:length(net.coeff), invcoeff{i} = inv(net.coeff{i}.'); end
    else
        invcoeff = {};
        coeff = net.coeff;
    end
    bvec = net.bvec;
    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        disp(['surrogate sample : ' num2str(k)]);
        S = single(Xorg);
        perm        = single(randperm(sigLen-lags));
        S(1:nodeNum,1:lags) = S(1:nodeNum,perm(1):perm(1)+lags-1); %Initialization of the AR surrogate
        perm2       = single(randperm(size(Err,2)));

        for t=lags+1:sigLen
            A = S(:,t);
            noise2 = noise(:,perm2(t-lags));
            for i=1:nodeNum
%            parfor i=1:nodeNum
                S2 = [];
                for p=1:lags
                    S2 = [S2; S(idxs{i,p},t-p)];
                end

                % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);
                % yield next time step
                if isempty(invcoeff)
                    score = (S2.' - mu{i}) / coeff{i}.';
                else
                    score = (S2.' - mu{i}) * invcoeff{i};
                end
                subScore = [score(:,1:mc{i}), 1];
                A(i) = subScore * bvec{i} + noise2(i);
            end
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

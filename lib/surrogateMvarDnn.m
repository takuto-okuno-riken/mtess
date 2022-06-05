%%
% Surrogate multivariate signal generation by multivariate VARDNN
% returns surrogated signals (Y)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mVARDNN network
%  dist         distribution of noise to yield surrogates ('gaussian'(default), 'residuals')
%  surrNum      output number of surrogate samples (default:1)
%  yRange       range of Y value (default:[Xmin-Xrange/5, Xmax+Xrange/5])

function Y = surrogateMvarDnn(X, exSignal, nodeControl, exControl, net, dist, surrNum, yRange)
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

    % calculate error residuals
    if isfield(net,'rvec')
        rvlen = length(net.rvec{1});
        for i=2:nodeNum, if rvlen > length(net.rvec{i}), rvlen = length(net.rvec{i}); end; end
        Err = single(nan(nodeNum,rvlen));
        for i=1:nodeNum
            Err(i,:) = net.rvec{i}(1:rvlen);
        end
    else
        Z = flipud(Xorg.'); % need to flip signal

        Err = nan(nodeNum,sigLen-lags);
        Zj = zeros(sigLen-lags, lags*inputNum);
        for p=1:lags
            Zj(:,1+inputNum*(p-1):inputNum*p) = Z(1+p:sigLen-lags+p,:);
        end
        nodeLayers = net.nodeLayers;
        nodeNetwork = net.nodeNetwork;
        for i=1:nodeNum
%        parfor i=1:nodeNum
            if isempty(nodeLayers), continue; end
            disp(['node residuals ' num2str(i)]);
            [~,idx] = find(control(i,:,:)==1);

            Xt = Z(1:sigLen-lags,i);
            Xti = Zj(:,idx);

            % predict and get errors
            A = predict(nodeNetwork{i}, Xti.', 'ExecutionEnvironment', 'cpu');
            Err(i,:) = Xt.' - A;
        end
        Z = [];  % clear memory
        Zj = []; % clear memory
    end

    idxs = {};
    for i=1:nodeNum
        for k=1:lags
            [~,idxs{i,k}] = find(control(i,:,k)==1);
        end
    end

    if strcmp(dist,'gaussian')
        P  = mean(Err.');
        EC = cov(Err.',1);
        noise = (mvnrnd(P,EC,size(Err,2)))';
    else
        noise = Err;
    end

    nodeNetwork = net.nodeNetwork;
    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        disp(['surrogate sample : ' num2str(k)]);
        S = Xorg;
        perm        = randperm(sigLen-lags);
        S(1:nodeNum,1:lags) = S(1:nodeNum,perm(1):perm(1)+lags-1); %Initialization of the AR surrogate
        perm2       = randperm(size(Err,2));

        for t=lags+1:sigLen
            A = S(:,t);
            noise2 = noise(:,perm2(t-lags));
            for i=1:nodeNum
%            parfor i=1:nodeNum
                if isempty(nodeNetwork{i}), A(i)=S(i,t-1); continue; end
                S2 = [];
                for p=1:lags
                    S2 = [S2; S(idxs{i,p},t-p)];
                end

                % yield next time step
                A(i) = predict(nodeNetwork{i}, S2, 'ExecutionEnvironment', 'cpu') + noise2(i);
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

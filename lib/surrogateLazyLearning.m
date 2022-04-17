%%
% Surrogate multivariate signal generation by Lazy Learning
% returns surrogated signals (Y)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  LL           LL structure
%  kn           number of k of nearest neighbor (default:1)
%  dist         distribution of noise to yield surrogates ('gaussian'(default), 'residuals')
%  surrNum      output number of surrogate samples (default:1)
%  yRange       range of Y value (default:[Xmin-Xrange/5, Xmax+Xrange/5])

function Y = surrogateLazyLearning(X, exSignal, nodeControl, exControl, LL, kn, dist, surrNum, yRange)
    if nargin < 9, yRange = NaN; end
    if nargin < 8, surrNum = 1; end
    if nargin < 7, dist = 'gaussian'; end
    if nargin < 6, kn = 1; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    inputNum = nodeNum + exNum;
    lags = LL.lags;
    if ~isempty(LL.Y)
        yLen = size(LL.Y,1);
    end

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
    for i=1:nodeNum
        for k=1:lags
            [~,idxs{i,k}] = find(control(i,:,k)==1);
        end
    end

    if ~isempty(LL.Y)
        Xt = {};
        for i=1:nodeNum
            Xti = [];
            for k=1:lags
                Xti = [Xti, LL.Y(1+k:yLen-lags+k,idxs{i,k})];
            end
            Xt{i} = Xti;
        end
        Yt = LL.Y(1:yLen-lags,:);
    else
        CY = LL.CY;
        cxNum = length(CY);
        idxn = cell(nodeNum,1);
        for n=1:nodeNum
            [~,idxn{n}] = find(control(n,:,:)==1);
        end
        allInLen = 0;
        for i=1:cxNum
            allInLen = allInLen + size(CY{i},1) - lags;
        end
        Xt = {};
        Xti = {};
        for n=1:nodeNum
%        parfor n=1:nodeNum
            Xt{n} = single(nan(allInLen,1));
            Xti{n} = single(nan(allInLen,length(idxn{n})));
            xts = 1;

            % this is redundant for every node. but it is necessary to avoid
            % too much memory consumption
            for i=1:cxNum
                % set node input
                Z = single(CY{i});

                sLen = size(Z,1);
                sl = sLen-lags;
                Yj = single(zeros(sl, lags*inputNum));
                for p=1:lags
                    Yj(:,1+inputNum*(p-1):inputNum*p) = Z(1+p:sl+p,:);
                end
                Xt{n}(xts:xts+sl-1,:) = Z(1:sl,n);
                Xti{n}(xts:xts+sl-1,:) = Yj(:,idxn{n});
                xts = xts + sl;
            end
            Z = [];  % clear memory
            Yj = []; % clear memory
        end
    end

    % calculate error residuals
    if ~isempty(LL.Y)
        Err = single(nan(nodeNum,sigLen-lags));
        for i=1:nodeNum
            % predict and get errors
            A = predictLazyLearning(Xt{i}.', Xt{i}.', Yt(:,i).', kn);
            Err(i,:) = Yt(:,i).' - A;
        end
    else
        if kn == 1
            Err = single(zeros(nodeNum,allInLen));
        else
            Err = single(nan(nodeNum,allInLen));
            for n=1:nodeNum
%            parfor n=1:nodeNum
                % predict and get errors
                A = predictLazyLearning(Xti{n}.', Xti{n}.', Xt{n}.', kn);
                Err(n,:) = Xt{n}.' - A;
            end
        end
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
            noise2 = noise(:,perm2(t-lags));
            if ~isempty(LL.Y)
                for i=1:nodeNum
%                parfor i=1:nodeNum
                    S2 = [];
                    for p=1:lags
                        S2 = [S2; S(idxs{i,p},t-p)];
                    end

                    % yield next time step
                    A(i) = predictLazyLearning(S2, Xt{i}.', Yt(:,i).', kn) + noise2(i);
                end
            else
                for n=1:nodeNum
            %    parfor n=1:nodeNum
                    S2 = [];
                    for p=1:lags
                        S2 = [S2; S(idxs{n,p},t-p)];
                    end

                    % yield next time step
                    A(n) = predictLazyLearning(S2, Xti{n}.', Xt{n}.', kn) + noise2(n);
                end
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

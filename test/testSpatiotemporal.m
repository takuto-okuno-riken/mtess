% test for Spatiotemporal model Shinn et al (2022)
% https://github.com/murraylab/spatiotemporal
% Compute TA-Δ1 (i.e. first-order temporal autocorrelation)
% Compute SA-λ and SA-∞ (i.e. measurements of spatial autocorrelation)

function testSpatiotemporal
    % check TA-Δ1
    TAs = [0, .2, .4, .6, .8];
    for i=1:length(TAs)
        load(['data/spatio/test_ta' num2str(TAs(i)) '.mat']);
        mta = calcTAD1(ts);
        if abs(mta - target_ta) >= .05
            disp(['error calcTAD1 ta=' num2str(mta) ', target_ta=' num2str(target_ta)]);
        end
        if abs(mta - ta) >= 1e-3
            disp(['error calcTAD1 ta=' num2str(mta) ', ta=' num2str(ta)]);
        end
    end

    % test SA-λ and SA-∞ 
    SAs = [10, 50, 5, 30];
    for i=1:length(SAs)
        load(['data/spatio/test_sa' num2str(SAs(i)) '.mat']);
        dists2 = pdist2(poss,poss);
        if any(abs(dists2(:) - dists(:)) >= 1e-3)
            disp(['error distance i=' num2str(i)]);
        end
        cm2 = spatialExponentialFloor(dists2, params(1), params(2));
        if any(abs(cm2(:) - cm(:)) >= 1e-3)
            disp(['error spatialExponentialFloor i=' num2str(i)]);
        end
        [SAlambda, SAinf] = calcSAlambda(cm2, dists2);
        if max(abs(log([SAlambda, SAinf]./params))) >= .1
            disp(['error calcSAlambda SA=' num2str(params(1)) ', SAlambda=' num2str(SAlambda)]);
        end
        if max(abs(log([SAlambda, SAinf]./fitparams))) >= .1
            disp(['error calcSAlambda SA=' num2str(SAs(i)) ', SAlambda=' num2str(SAlambda) ', target_SAlambda=' num2str(fitparams(1))]);
        end
    end

    % test eigen matrix surrogate
    for i=1:2
        load(['data/spatio/test_egmsurr' num2str(i-1) '.mat']);
        cm2 = calcFunctionalConnectivity(ts);
        if any(abs(cm2(:) - cm(:)) >= 1e-3)
            disp(['error FC i=' num2str(i)]);
        end
        e2 = flip(svd(cm2)); % eig does not match scipy.linalg.eigvalsh(cm), we use svd instead.
        if any(abs(ev_cm(:) - e2(:)) >= 1e-3)
            disp(['found diff between eig i=' num2str(i)]);
        end
        Y = surrogateEigenMatrix(cm2);
        e3 = flip(svd(Y));
        if max(abs(log(e2./e3))) >= 1e-3
%        if any(abs(e3(:) - e2(:)) >= .01)
            disp(['error surrogateEigenMatrix i=' num2str(i)]);
        end
    end

    % test eigen time-series surrogate
    for i=1:2
        load(['data/spatio/test_egtsurr' num2str(i-1) '.mat']);
        cm2 = calcFunctionalConnectivity(ts);
        if any(abs(cm2(:) - cm(:)) >= 1e-3)
            disp(['error FC i=' num2str(i)]);
        end
        e2 = flip(svd(cm2));
        if any(abs(ev_cm(:) - e2(:)) >= 1e-3)
            disp(['found diff between eig i=' num2str(i)]);
        end
        X = surrogateEigenTS(ts);
        cm3 = calcFunctionalConnectivity(X);
        e3 = flip(svd(cm3));
        if mean(abs(log(e2./e3))) >= .5
            disp(['error surrogateEigenTS i=' num2str(i), ', max_err=' num2str(max(abs(log(e2./e3))))]);
        end
    end

    % test spaciotemporal surrogate
    tsLen = 1000; % Big to get a better fit
    tr = 1; % 1 sec
    for i=1:3
        load(['data/spatio/test_sptmsurrno' num2str(i-1) '.mat']);
        params = double(params); % hmm.
        % surrogate spatio temporal time-series
        X = surrogateSpatiotemporal(distance_matrix, params(1), params(2), ta_delta1s, tsLen, tr, 1);
        figure; plot(X'); title(['surrogate spatio temporal saLambda=' num2str(params(1)) ', saInf=' num2str(params(2))]);

        % check TA-Δ1 with surrogate time-series
        mta = calcTAD1(X);
        if max(abs(newtas-mta')) >= .3 % check with python implementation
            disp(['error surrogateSpatiotemporal TA-Δ1 with python, i=' num2str(i), ', max_err=' num2str(max(abs(newtas-mta')))]);
        end
        if max(abs(ta_delta1s-mta')) >= .3 % check with original
            disp(['error surrogateSpatiotemporal TA-Δ1 with original, i=' num2str(i), ', max_err=' num2str(max(abs(ta_delta1s-mta')))]);
        end

        % check SA-λ and SA-∞ with surrogate time-series
        cm3 = calcFunctionalConnectivity(X);
        [SAlambda, SAinf] = calcSAlambda(cm3, distance_matrix);
        if params(2)==0, params(2)=nan; SAinf=nan; end % ignore params(2)==0 case.
        if max(abs(log(fitparams./[SAlambda, SAinf]))) >= .3 % check with python implementation
            disp(['error surrogateSpatiotemporal SA with python, SA=' num2str(params(1)) ', SAlambda=' num2str(SAlambda) ', SAinf=' num2str(SAinf)]);
        end
        if max(abs(log(params./[SAlambda, SAinf]))) >= .3 % check with original. 
            disp(['error surrogateSpatiotemporal SA with original, SA=' num2str(params(1)) ', SAlambda=' num2str(SAlambda)]);
        end
    end
end

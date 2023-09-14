%%
% Surrogate spatiotemporal time-series (Shinn et al., 2022)
% returns surrogated signals Y (node x time seriese x surrNum)
% input:
%  D           distance matrix (node x node)
%  saLambda    SA-λ value of spatial autocorrelation
%  saInf       SA-∞ value of spatial autocorrelation
%  TA          TA-Δ1 (first-order temporal autocorrelation) (node x 1 values)
%  tsLen       length of timeseries to generate
%  tr          the spacing between timepoints (sec)
%  highpass    highpass filer (Hz) (not supported)
%  surrNum     output number of surrogate samples (default:1)

function Y = surrogateSpatiotemporal(D, saLambda, saInf, TA, tsLen, tr, highpass, surrNum)
    if nargin < 8, surrNum = 1; end
    if nargin < 7, highpass = 0; end
    nodeNum = size(D,1);

    % Determine the pink noise exponent alpha from the TA-delta1
    alphas = zeros(1, length(TA));
    for i=1:length(TA)
        alphas(i) = TA2alpha(tsLen, tr, max(0,TA(i)), highpass);
    end
    % Use these alpha values to construct desired frequency spectra
    spectra = zeros(length(alphas),floor(tsLen/2)+1);
    for i=1:length(alphas)
        spectra(i,:) = makeSpectrum(tsLen, tr, alphas(i), highpass);
    end
    % Spatial autocorrelation matrix
    cm = spatialExponentialFloor(D, saLambda, saInf);
    % Compute timeseries from desired correlation matrix and frequency spectra
    Y = nan(nodeNum,tsLen,surrNum);
    for i=1:surrNum
        Y(:,:,i) = spatialTemporalTimeseries(cm, spectra);
    end
end

function TS = spatialTemporalTimeseries(CM, spectra)
    nodeNum = size(CM,1);
    nFreqs = size(spectra,2);
    tsLen = (nFreqs-1)*2;
    SS = sum(spectra.^2, 2);
    CS = (spectra * spectra') ./ sqrt(SS * SS'); % cosine similarity
    CV = CM ./ CS; % covariance matrix
    % must be a positive semi-definite matrix.
    [V,D] = eig(CV);
    E = diag(D);
    if min(E) < 0
        E(E < 0) = 0;
        s = sum(E);
        E(end) = E(end) - (s - length(E));
        disp(['Warning (spatialTemporalTimeseries): eigenvalues of covmat were less than zero in source matrix by ' num2str(s - length(E))]);
        CV = V * diag(E) * inv(V);
    end
    M  = zeros(1,nodeNum);
    rvs = mvnrnd(M,CV,nFreqs*2);
    reals = rvs(1:nFreqs,:)' .* spectra;
    ims = rvs(nFreqs+1:end,:)' .* spectra;
    % Since the signal length is even, frequencies +/- 0.5 are equal
    % so the coefficient must be real.
    ims(:,end) = 0;
    % The DC component must be zero and real
    reals(:,1) = 0;
    ims(:,1) = 0;
    TS = real(ifft(reals + 1i*ims,tsLen,2));
end

function spectrum = makeSpectrum(tsLen, tr, alpha, highpass)
    if nargin < 4, highpass = 0; end

%   freqs = np.fft.rfftfreq(tslen, tr) in Python
    val = 1.0/(tsLen*tr);
    N = floor(tsLen/2) + 1;
    freqs = 0:1:N-1;
    freqs = freqs * val;
    spectrum = freqs.^(-alpha/2);

    % highpass mode is not supported
    if highpass > 0
    end

    spectrum(1) = 0;
end

function ta = getSpectrumTA(spectrum)
    n = length(spectrum);
    weightedsum = sum(spectrum.^2 .* cos(pi * (0:1:n-1)/n));
    ta = weightedsum / sum(spectrum.^2);
end

function x = TA2alpha(tslen, tr, target_ta, highpass)
    v0 = [1.5]; % initial value
    fun = @(v)(getSpectrumTA(makeSpectrum(tslen, tr, v(1), highpass)) - target_ta)^2; % error evaluation function
%    options = optimset('Display','iter');
    [x,fval] = fminsearch(fun,v0); %,options);
end

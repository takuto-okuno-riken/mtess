%%
% Get auto range for surrogate signal
% returns [min max] (yRange)
% input:
%  X            multivariate time series matrix (node x time series)

function yRange = getAutoRange(X)
    xmax = max(X,[],'all');
    xmin = min(X,[],'all');
    width = xmax - xmin;
    yRange = [xmin - width/4, xmax + width/4];
end

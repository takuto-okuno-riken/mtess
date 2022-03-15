%%
% Plot MTESS all matrix
% input:
%  P                MTESS statistical properties matrix (cell number x 7)
%  ccLags           time lags for Cross-Correlation function (default: 8)
%  pccLags          time lags for Partial Cross-Correlation function (default: 8)
%  pcName           Partial Correlation function name (default: 'PC')
function plotMtessSpiderPlot(P, ccLags, pccLags, pcName)
    if nargin < 4, pcName = 'PC'; end
    if nargin < 3, pccLags = 8; end
    if nargin < 2, ccLags = 8; end
    spider_plot(P,'AxesLimits',[0,0,0,0,0,0,0;5,5,5,5,5,5,5],'AxesPrecision',[0], 'AxesDisplay', 'one', ...
        'AxesLabels', {'mean','SD','Amp','FC',pcName,['CC' num2str(ccLags)],[pcName 'C' num2str(pccLags)]});
end

%%
% Plot MTESS all matrix
% input:
%  P                MTESS statistical properties matrix (cell number x 8)
function plotMtessSpiderPlot(P)
    spider_plot(P,'AxesLimits',[0,0,0,0,0,0,0,0;5,5,5,5,5,5,5,5],'AxesPrecision',[0], 'AxesDisplay', 'one', ...
        'AxesLabels', {'SD','AC','PAC','CM','PCM','CCM','PCCM','mKT'});
end

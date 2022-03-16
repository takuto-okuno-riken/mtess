%%
% Plot MTESS all matrix
% input:
%  MTS              MTESS matrix (cell number x cell number)
%  MTSp             MTESS statistical properties matrix (cell number x cell number x 7)

function plotMtessAllMatrix(MTS, MTSp, outname)
    if nargin < 3, outname = ''; end

    % subject vs. subject all matrix
    mt = MTS; sc = MTSp;
    mt(isnan(mt)) = 0;
    sc(isnan(sc)) = 0;
    S(:,:,1) = sc(:,:,3).' + mt;
    S(:,:,2) = sc(:,:,2).' + sc(:,:,1);
    S(:,:,3) = sc(:,:,5).' + sc(:,:,4);
    S(:,:,4) = sc(:,:,7).' + sc(:,:,6);
    
    outname = strrep(outname,'_','-'); % for text bug.
    prop2 = {'MTESS/Amp','Mean/SD','FC/PC','CC/PCC'};
    for i=1:4
        % show MTESS matrix
        figure;
        clims = [0,5];
        imagesc(S(:,:,i),clims);
        daspect([1 1 1]);
        title([prop2{i} ' matrix of ' outname]);
        xlabel('Cell number');
        ylabel('Cell number');
        colorbar;
    end
end

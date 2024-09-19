%%
% Adjust NIfTI volume direction based on NIfTI info Transpose matrix
% returns adjusted NIfTI volume
% input:
%  V            nifti 4D volume (X x Y x Z x frames)

function V = adjustVolumeDir(V, info)
    if strcmp(info.TransformName,'Qform')
        disp(['NIfTI Qform format is not supported : ' info.Filename])
        return;
    end
    A = [0 1 1; 1 0 1; 1 1 0];
    if sum(abs(A .* info.Transform.T(1:3,1:3)),'all') > 0
        disp(['This transformation is not supported : ' info.Filename])
        return;
    end
    if info.Transform.T(1,1) < 0   % check flip X axis
        V = flipud(V);
    end
    if info.Transform.T(2,2) < 0   % check flip Y axis
        V = fliplr(V);
    end
    if info.Transform.T(3,3) < 0   % check flip Z axis
        V = flip(V,3);
    end
end


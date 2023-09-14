%%
% Surrogate matrix by eigenvalue (compatible with scipy.stats.random_correlation in Python) (Shinn et al., 2022)
% returns surrogated matrix Y (node x node x surrNum)
% input:
%  C          correlation matrix (node x node)
%  surrNum    output number of surrogate samples (default:1)

function Y = surrogateEigenMatrix(C, surrNum)
    if nargin < 2, surrNum = 1; end
    n = size(C,1);
    [V,D] = eig(C);
    E = diag(D);
    if min(E) < 0
        E(E < 0) = 0;
        s = sum(E);
        E(end) = E(end) - (s - length(E));
        disp(['Warning: eigenvalues were less than zero in source matrix by ' num2str(s - length(E))]);
    end
    D = diag(E);
    Y = nan(n,n,surrNum);
    for i=1:surrNum
        p = randperm(n);
        S = V(p,p); % permutate original matrix (not fully random)
        Yi = S*D*inv(S);
        Yi = to_corr_(Yi); % from SciPy
        Y(:,:,i) = Yi;
    end
end

% scipy.stats.random_correlation in Python can get diagonal==1 correlation matrix.
% But it requires careful rotation by scipy.stats._multivariate.random_correlation_gen._to_corr()
function M = to_corr_(M)
    d = size(M,1);
    for i=1:d-1
        if M(i,i) == 1
            continue;
        elseif M(i,i) > 1
            for j=i+1:d
                if M(j,j) < 1, break; end
            end
        else
            for j=i+1:d
                if M(j,j) > 1, break; end
            end
        end

        [c, s] = givens_to_1_(M(i,i), M(j,j), M(i,j));

        G = eye(d);
        G(i,i) = c; G(j,j) = c;
        G(j,i) = -s; G(i,j) = s;
        M = G' * M * G;
    end
end

function [c, s] = givens_to_1_(aii, ajj, aij)
    aiid = aii - 1;
    ajjd = ajj - 1;
    if ajjd == 0
        c=0; s=1;
        return;
    end

    dd = sqrt(max([aij^2 - aiid*ajjd, 0]));
    t = (aij + dd * sign(aij)) / ajjd;
    c = 1 / sqrt(1 + t*t);
    if c == 0
        s = 1.0; % Underflow
    else
        s = c*t;
    end
end

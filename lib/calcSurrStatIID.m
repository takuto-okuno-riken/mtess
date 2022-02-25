%%
% discriminating statistic function to check IID (independent and identically distributed) for surrogate rank test
% returns discriminating statistic value (t)
% based on J.Theilear and D.Prichard, Physica D: Nonlinear Phenomena (1996) pp.221-235.
% input:
%  X          time series vector
%  params     discriminating statistic function parameters

function t = calcSurrStatIID(X, params)
    n = length(X);
    m = mean(X);
    X0 = X(1:end-1) - m;
    X1 = X(2:end) - m;
    X2 = X0(:).' * X1(:) / (n - 1);
    m2 = mean((X - m).*(X - m));
    t = X2 / m2;
end

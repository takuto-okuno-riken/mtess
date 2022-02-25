%%
% discriminating statistic function to check linearity for surrogate rank test
% returns discriminating statistic value (t)
% based on J.Theilear and D.Prichard, Physica D: Nonlinear Phenomena (1996) pp.221-235.
% input:
%  X          time series vector
%  params     discriminating statistic function parameters

function t = calcSurrStatLinear(X, params)
    n = length(X);
    m = mean(X);
    m2 = mean((X - m).*(X - m));
    t1 = nan(1,n-2);
    for i=3:n
        t1(i-2) = X(i) - f(X(i-1), X(i-2), params);
    end
    t2 = mean(t1 .* t1);
    t = t2 / m2;
end

function w = f(u, v, c)
    w = c{1} + c{2}*u + c{3}*v + c{4}*u*u + c{5}*u*v + c{6}*v*v;
end

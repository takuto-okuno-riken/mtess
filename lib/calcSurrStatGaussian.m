%%
% discriminating statistic function to check gaussian for surrogate rank test
% returns discriminating statistic value (t)
% based on J.Theilear and D.Prichard, Physica D: Nonlinear Phenomena (1996) pp.221-235.
% input:
%  X          time series vector
%  params     discriminating statistic function parameters

function t = calcSurrStatGaussian(X, params)
    m = mean(X);
    x2 = (X - m);
    m2 = mean(x2 .* x2 .* x2 .* x2);
    m3 = mean(x2 .* x2);
    t = m2 / (m3 * m3);
end

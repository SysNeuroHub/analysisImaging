function [a,b,yfit] = fitExp1D(y)
% fitExp1D fits an exponential curve y = a*exp(b*x) to a 1D vector
%
% INPUT
% y : 1D vector (row or column)
%
% OUTPUT
% a, b : fitted parameters
% yfit : fitted curve values

    y = y(:);                 % ensure column vector
    x = (1:length(y))';       % index as x-axis

    % remove non-positive values (log requires positive)
    idx = y > 0;
    x = x(idx);
    y = y(idx);

    % linearize: log(y) = log(a) + b*x
    p = polyfit(x, log(y), 1);

    b = p(1);
    a = exp(p(2));

    % compute fitted curve for full x
    xfull = (1:length(y))';
    yfit = a * exp(b * xfull);
end
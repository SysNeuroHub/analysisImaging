function [k,beta, meany] = fitStretchedExp(y,t,estimateBeta)

opts = optimoptions('lsqnonlin',...
    'Display','off',...
    'MaxIterations',200);


if estimateBeta

    p0 = [5e-4 0.7];

    fun = @(p) exp(-p(1)*(t.^p(2))) - y;

    p = lsqnonlin(fun,...
        p0,...
        [0 0.2],...
        [0.05 1],...
        opts);

    k = p(1);
    beta = p(2);


else

    beta = estimateBeta;

    fun = @(kk) exp(-kk*(t.^beta))-y;

    k = lsqnonlin(fun,...
        5e-4,...
        0,...
        0.05,...
        opts);

end

end
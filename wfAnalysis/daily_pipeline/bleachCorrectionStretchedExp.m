function [Icorr,param,diag] = bleachCorrectionStretchedExp(I,betaFixed)

if nargin < 1
    error('Input matrix I is required.');
end


I = double(I);

[nTrace,nT] = size(I);

tFit = (1:nT)';


if nargin == 1

    estimateBeta = true;
    beta = zeros(nTrace,1);

else

    estimateBeta = false;
    beta = repmat(betaFixed,nTrace,1);

end

k = zeros(nTrace,1);

for ii=1:nTrace %parfor not working

    y = I(ii,:)';

    scale = median(y(1:min(20,nT)));

    if scale<=0 || ~isfinite(scale)
        continue
    end

    y = y/scale;


    if estimateBeta
        [k(ii),beta(ii)] = fitStretchedExp(y,tFit,true);
    else
        [k(ii),~] = fitStretchedExp(y,tFit,betaFixed);
    end

end


param.k = k;
param.beta = beta;
param.scaleFac = [];

%% Generate correction without scaling
Icorr_unscaled = applyStretchedExp(I,param);

% Preserve global mean
param.scaleFac = mean(I(:),'omitnan') / mean(Icorr_unscaled(:),'omitnan');

%% Generate correction with scaling
Icorr = applyStretchedExp(I,param);


diag.meanBefore = mean(I(:),'omitnan');
diag.meanAfter = mean(Icorr(:),'omitnan');
diag.meanSignal = mean(Icorr,1);


end



%% =========================================================
% Local functions
%% =========================================================

function r = stretchedResidual(p,t,y)

r = exp(-p(1)*(t.^p(2))) - y;

end


function r = stretchedResidualFixed(k,t,beta,y)

r = exp(-k*(t.^beta)) - y;

end
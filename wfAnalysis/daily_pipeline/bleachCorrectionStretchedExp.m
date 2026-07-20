function [Icorr,param,diag] = bleachCorrectionStretchedExp(I,betaFixed,frameNumbers)

if nargin < 1
    error('Input matrix I is required.');
end


I = double(I);

[nTrace,nT] = size(I);

if nargin < 3
    frameNumbers = (1:nT)';
end

% Progress bar
N = nTrace;
count = 0;


% D = parallel.pool.DataQueue;
% afterEach(D,@updateWaitbar);


if nargin == 1
    
    estimateBeta = true;
    beta = zeros(nTrace,1);
    
else
    
    estimateBeta = false;
    beta = repmat(betaFixed,nTrace,1);
    
end

k = zeros(nTrace,1);

if estimateBeta
    parfor ii=1:nTrace
        
        y = I(ii,:)';
        
        scale = median(y(1:min(20,nT)));
        
        if scale<=0 || ~isfinite(scale)
            continue
        end
        
        y = y/scale;
       
        [k(ii),beta(ii)] = fitStretchedExp(y,frameNumbers,true);
        
         % Notify completion of one iteration
         % send(D,1);
    end
    
else %% estimate k but not beta
    parfor ii=1:nTrace
        y = I(ii,:)';
        
        scale = median(y(1:min(20,nT)));
        
        if scale<=0 || ~isfinite(scale)
            continue
        end
        
        y = y/scale;
        
        [k(ii),~] = fitStretchedExp(y,frameNumbers,betaFixed);
    
        % Notify completion of one iteration
        % send(D,1);
    end
end

% close(h);

param.k = k;
param.beta = beta;
param.scaleFac = [];

%% Generate correction without scaling
Icorr_unscaled = applyStretchedExp(I,param);

% Preserve global mean

%% Generate correction with scaling
param.scaleFac = mean(I,2,'omitnan') ./ mean(Icorr_unscaled,2,'omitnan');
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

function updateWaitbar(~)
    count = count + 1;
    if mod(count,100)==0 || count==N
        fprintf('%d/%d (%.1f%%)\n',count,N,100*count/N);
    end
end
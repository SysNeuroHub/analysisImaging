function [ampMap, map, combMap] = retinotopyMakeMap(rrr,p)
%
% [ampMap map combMap] = retinotopyMakeMap(rrr,p)
%
% rrr (nr nc istim), p is from ProtocolLoad
%
% pass -rrr for Intrinsic (activations <0)
%
%
% the logic of the following is:
% given data r_i and model k*m_i, the mean squared error is
% sum_i((k*m_i - r_i)^2) = sum_i(ri^2) + k^2*sum_i(mi^2) -2*k*sum_i(m_i*r_i)
% which in vector form is RR + k^2 *MM - 2*k*MR
%

[nr, nc, nstim] = size(rrr); 
first_stim = 1; 
last_stim = nstim;
Selection = 2;
np = nr*nc;

rr = reshape(rrr, np, nstim);

% the values of sigma that we will try
sigmas = 0.5:0.1:3.5;

nsigmas = length(sigmas);

mm = cell(nsigmas,1);

iii = cell(nsigmas,1);
kkk = cell(nsigmas,1);
MedianErrs = nan(nsigmas,1);

figure; % to see the plot updating as the computation goes on

for isigma = 1:nsigmas
    
    sigma = sigmas(isigma); % should optimize here
    
    mm{isigma} = zeros(nstim,50);
    PrefPositions = linspace(-0.5,nstim+1.5,50);
    for ipref = 1:50
        mm{isigma}(:,ipref) = gaussian( [PrefPositions(ipref), 1, 0, sigma], 1:nstim );
    end
    
    RM = rr*mm{isigma};                         % np X 50
    RR = repmat(sum(rr.^2,2),[1 50]);   % np X 50
    MM = repmat(sum(mm{isigma}.^2,1),[np 1]);   % np X 50
    
    % for each of the np and each of the 50, find the best k
    
    K = RM ./ MM; % np X 50
    
    K = max(K,0);
    
    E = K.^2 .* MM + RR - 2* K .*RM;
    
    [ee, ii] = min(E,[],2);
    eee = reshape(ee, [nr,nc]);
    
    MedianErrs(isigma) = median(eee(any(rrr>prctile(rrr(:),90),3)));
    
    %     figure; imagesc( eee ); axis image; colorbar; colormap bone;
    %     title(sprintf('Error with sigma = %2.2f is %2.1d',sigma,MedianErrs(isigma)));
    %     drawnow
    
    fprintf('Sigma = %2.2f (%d of %d) deviation %2.1d\n', sigma, isigma, nsigmas,MedianErrs(isigma));
    
    kk = zeros(np,1);
    for ip = 1:np
        kk(ip) = K(ip,ii(ip));
    end
    
    iii{isigma} = reshape(ii,[nr,nc]);
    kkk{isigma} = reshape(kk,[nr,nc]);
    
    plot(sigmas,MedianErrs); 
    xlabel('Sigma');
    ylabel('Median deviation');
    drawnow;
end

[BestErr, iBest] = min(MedianErrs);
if iBest == 1 || iBest==nsigmas
    disp('Could not find best sigma');
    sigma = 1;
else
    sigma = sigmas(iBest);
end

hold on
plot( sigmas([1 end]), BestErr*[1 1], 'k--' );

figure('Name','Tuning Strength');
imagesc( kkk{iBest} ); axis image; colorbar; colormap bone
title('Map of tuning strength');

map = PrefPositions(iii{iBest});
try
    PrefPositionRealValues = linspace(p.pars(p.activepars{Selection},1),p.pars(p.activepars{Selection},nstim),50);
catch
    PrefPositionRealValues = linspace(p.pars(p.activepars{1}(end),1),p.pars(p.activepars{1}(end),nstim),50); %May 11 AP for experiments older than Oct 2010
end
mapRealValues = PrefPositionRealValues(iii{iBest});

figure('Name','Retinotopy');
imagesc(map); axis image; 
cb = colorbar;
try
    StimValues = p.pars(p.activepars{Selection},first_stim:last_stim); %AP Apr 11 to deal with multiple activepar
catch
    StimValues = p.pars(p.activepars{1}(end),:);%
    % first_stim, last_stim are newly introduced by AP?
end
try
    ylabel(cb,p.pardefs{p.activepars{Selection}});
catch% some experiments can not be loaded here. 20110422 TS
    %for experiments older than Oct 2010 use instead --- May 2011 AP 
    ylabel(cb,p.pardefs{p.activepars{1}(end)});
end
set(cb,'ytick',1:nstim,'yticklabel',StimValues);
title('Map of preferred positions');

%% Combine the two

% combine position and amplitude maps
ampMap  = kkk{iBest};
retMap  = map;

nampMap = nthroot(ampMap,4);
nampMap = nampMap - min(nampMap(:));
nampMap = nampMap./max(nampMap(:));

rgbMap = ind2rgb( gray2ind( mat2gray( retMap ), 64 ), jet);
combMap = []; 
for ii = 1:3
    combMap(:,:,ii) = squeeze(rgbMap(:,:,ii)).*nampMap(:,:);
end

% plot
figure('Name','Retinotopy'); clf
imagesc(combMap); axis image; % axis xy
cb = colorbar;

try
    StimValues = p.pars(p.activepars{Selection},first_stim:last_stim); %AP Apr 11 to deal with multiple activepar
catch
    StimValues = p.pars(p.activepars{1}(end),:);%
    % first_stim, last_stim are newly introduced by AP?
end

ylabel(cb,p.pardefs{p.activepars{1}(end)});
set(cb,'ytick',1:nstim,'yticklabel',StimValues);
title('Map of preferred positions');



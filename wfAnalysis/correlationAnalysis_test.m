%
%TODO
% global signal regression?
% automatic cell detection via s2p?
% estimate coupling via GLM Runyan ... Harvey 2017 Harvey


addpath(genpath('C:\npy-matlab'));
addpath(genpath('C:\Users\dshi0006\npy-matlab'));
addpath(genpath('C:\Users\Analysis\npy-matlab'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\widefield'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\visbox\rigbox'));
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis');
addpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis');


%% experiment
%ppbox notation
expt.subject = 'TIGRE2GCaMP6s_318';
expt.expDate = '2020-11-09_1';
expt.expNum = 1;
bklightCtrl = 0;


%% SVD
nSV = 1000;
params.movieSuffix = 'blue';% 'purple''corr_dFF'; %'blue'  %U for corr_dFF can have NANs..
params.useCorrected = 1;

%% analysis
%highpassCutoff = 0.01; %[Hz]
%lowpassCutoff = []; %[Hz]
resizeS = 1;%0.25; %spatial rescaling factor
% pixies = [169 174; 308 197; 263 150; 127 136; 149 252];
% pixies = [285 167; ...
%     295 139; ...
%     95 150; ...
%     55 137; ...
%     113 162; ...
%     92 186; ...
%     188 144; ...
%     133 179; ...
%     216 191; ...
%     155 169; ...
%     124 260;...
%     188 244;...
%     105 116;...
%     194 235];
% pixies = 2*[129 195;...
%     149 290; 
%     155 410; ...
%     342 92;...
%     259 555; ...
%     406 312; ...
%     310 443;  ...
%     360 385];
pixies = [31 282;...
    47 222;...
    103 284;...
    112 242;...
    128 189;...
    139 262;...
    215 252;...
    259 234;...
    300 322;...
    411 271;...
    427 225];
    
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
figname = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
    params.movieSuffix];



resultSaveDir = fullfile('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\AnalysisResult\wf',...
    [expt.subject '_' expt.expDate '_' num2str(expt.expNum)]);
mkdir(resultSaveDir);


%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;
mpepDir = dat.reposPath('main', 'master');


%% load wf data
disp('Loading widefield data');
disp(expt)
[U, V, t, mimg, mask] = quickLoadUVt(expPath, nSV, saveVpath, params);
Fs = 1/median(diff(t));

%% long traces
Uflat = reshape(U, size(U,1)*size(U,2), size(U,3));
pixies_ind = sub2ind([size(U,1) size(U,2)], pixies(:,1), pixies(:,2));
ROItrace = (Uflat(pixies_ind,:)*V)';
stackedplot(t, ROItrace,'displaylabels',num2cell(1:size(pixies,1)),...
    'gridvisible','on');
 saveas(gcf, fullfile(resultSaveDir, [figname 'longTraces' ...
        num2str(size(pixies,1)) '.png']));
close(gcf);

%% seed-pixel correlation map
for ff = 1:size(pixies,1)
    pixelCorrelationViewerSVD(U, V,[pixies(ff,1) pixies(ff,2)]);
    set(gcf, 'position',[0 0 400 400]);
    
    saveas(gcf, fullfile(resultSaveDir, [figname 'seedPixCorr_at' ...
        num2str([pixies(ff,1) pixies(ff,2)]) '.png']));
    close(gcf);
end


%% compute seed-pixel correlation (taken from pixelCorrelationViewerSVD)
varCalcMax = false;
Ur = reshape(U, size(U,1)*size(U,2),[]); % P x S
covV = cov(V'); % S x S % this is the only one that takes some time really
varP = dot((Ur*covV)', Ur'); % 1 x P
ySize = size(U,1); xSize = size(U,2);

corrMat_s = zeros(size(pixies,1));
for ff = 1:size(pixies,1)
    pixelInd = sub2ind([ySize, xSize], pixies(ff,1), pixies(ff,2));
    covP = Ur(pixelInd,:)*covV*Ur'; % 1 x P
    if varCalcMax
        stdPxPy = varP(pixelInd).^0.5 * max(varP(:)).^0.5; % 1 x P
    else
        stdPxPy = varP(pixelInd).^0.5 * varP.^0.5; % 1 x P
    end
    corrMat = covP./stdPxPy; % 1 x P
    corrMat = reshape(corrMat, ySize, xSize);
    
    for ff_tgt = 1:size(pixies,1)
        corrMat_s(ff,ff_tgt) = corrMat(pixies(ff_tgt,1), pixies(ff_tgt,2));
    end
end

figure;
%% correlation matrix between selected pixels
clim = [0 1];
cmap = 'jet';

ax1=subplot(121);
imagesc(corrMat_s);axis equal tight;
set(ax1,'xtick',1:size(pixies,1),'ytick',1:size(pixies,1));
caxis(clim);
colormap(ax1,cmap);
mcolorbar(ax1,.5);

%% show the correlation matrix on the cortex
ax2=subplot(122);
imagesc(mimg); colormap(ax2,'gray');axis equal tight

%plot(pixies(:,1), pixies(:,2), 'ro');
hold on;
matcolors = value2Color(corrMat_s,clim,cmap);
for ff = 1:size(pixies,1)
    for ff_tgt = 1:size(pixies,1)
        line([pixies(ff,2) pixies(ff_tgt,2)],[pixies(ff,1) pixies(ff_tgt,1)], ...
            'color', matcolors(ff,ff_tgt,:));
            %             'linewidth',abs(corrMat_s(ff,ff_tgt)));
            hold on
    end
end
text(pixies(:,2), pixies(:,1), num2cell(1:size(pixies,1)), 'color', 'w');
saveas(gcf, fullfile(resultSaveDir, [figname 'seedPixMatrix_' ...
    num2str(size(pixies,1)) 'seeds.png']));
close(gcf);
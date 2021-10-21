function hemoCorrectSequential(Exps, pixSpace, FreqRange, vidNames_regressee, ...
    vidNames_regressor, vidName_align, regressFiltered)
% hemoCorrectSequential(Exps, pixSpace, FreqRange, vidNames_regressee, ...
%    vidNames_regressor, vidName_align, regressFiltered)

% created from hemoCorrectSequential for Shimaoka 2019 elif

if nargin < 7
    regressFiltered = 0;
end

nRegressors = length(vidNames_regressor);
if size(FreqRange,1) == 1
    for ch2 = 1:nRegressors
        FreqRange(ch2,:) = FreqRange(1,:);
    end
end
scrsz = get(0,'screensize');

disp('hemoCorrectSequential: loading U...');
thisDate = ['20' num2str(Exps.iseries(3:4)) '-' num2str(Exps.iseries(5:6)) ...
    '-' num2str(Exps.iseries(7:8))];
thisPath = fullfile('\\zserver\Data\Subjects', Exps.animal, thisDate);
U_common = readNPY([thisPath '\svdSpatialComponents' vidName_align '.npy']);

load( [thisPath '\dataSummary' vidName_align '.mat'] );
Fs = 1/median(diff(dataSummary.timeStampsFromStamp));

    
%% load dV of regressors
disp('hemoCorrectSequential: loading regressors...');
for ch2 = 1:length(vidNames_regressor)
    dVname = ['svdTemporalComponents' vidNames_regressor{ch2} '_BY' vidName_align];
    try
        dV = readNPY([thisPath '\' num2str(Exps.iexp) '\' dVname '.npy']);
    catch err
        dV = readNPY([thisPath '\' dVname '.npy']);
    end
    dVaux{ch2} = dV';
    clear dV
end

dVout2 = [];
for ch = 1:length(vidNames_regressee)
    %% load dV of regressee
    dVname = ['svdTemporalComponents' vidNames_regressee{ch} '_BY' vidName_align];
    try
        dVin = readNPY([thisPath '\' num2str(Exps.iexp) '\' ...
            dVname '.npy']);
        Vpath = [thisPath '\' num2str(Exps.iexp)];
    catch err
        dVin = readNPY([thisPath '\' dVname '.npy']);
        Vpath = thisPath;
    end
    dVin = dVin';
    
    %% sequential regression of dVin by dVaux{1} then dVaux{2}
    disp('hemoCorrectSequential: computing regression...');
    for ch2 = 1:nRegressors
        if ch2 > 1
            dVin = dVout{ch2-1};
        end
        [dVout{ch2},sfac{ch2}] = HemoCorrectLocal_simple(U_common, dVin, dVaux{ch2}, ...
            Fs, FreqRange(ch2,:), pixSpace, regressFiltered);
    end
    
    %[dVout2,sfac2] = HemoCorrectLocal_simple(U_common, dVout, dVaux{2}, Fs, FreqRange, pixSpace, regressFiltered);
    
    %     vidName_target = [vidNames_regressee{ch} '_reg' vidNames_regressor{1} ...
    %         vidNames_regressor{2}];
    vidName_target = [vidNames_regressee{ch} '_regressed'];
    if regressFiltered 
        vidName_target = [vidName_target '_f'];
    end
        disp(['hemoCorrectSequential: saving ' vidName_target]);
    saveV(dVout{end}, [], Vpath, vidName_target, vidName_align);
    
    %% plot scale factors
    h=[];
    for ch2 = 1:length(vidNames_regressor)
        sfac{ch2}(isinf(sfac{ch2}))=0;
        
        subplot(1,length(vidNames_regressor), ch2);
        imagesc(sfac{ch2});axis equal tight
        title(['regression coef by ' vidNames_regressor{ch2}]);
        caxis(prctile(sfac{ch2}(:),[1 99]));
        colorbar;
    end       
    set(gcf,'position',scrsz);
    screen2png(fullfile(Vpath, ['scalefactors_in' vidName_target]));
    close all
    clear dVin
    
    save(fullfile(Vpath, ['svdTemporalComponents' vidName_target '_BY' vidName_align]),...
        'sfac','FreqRange','pixSpace');
    
end

function saveV(svdTemporalComponents, t, Vpath, vidName_target, vidName_align)
%taken from saveSVD

fn = fullfile(Vpath, ['svdTemporalComponents' vidName_target '_BY' vidName_align]);
fnT = fullfile(Vpath, ['svdTemporalComponents' vidName_target '_BY' vidName_align '.timestamps']);

% if isfield(ops, 'saveAsNPY') && ops.saveAsNPY
writeUVtoNPY([], svdTemporalComponents, [], fn);
if ~isempty(t)
    writeNPY(t, [fnT '.npy']);
end
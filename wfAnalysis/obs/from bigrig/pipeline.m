addpath('\\zouni\Users\Experiment\Documents\MATLAB\dsbox\SVDtools\'); %VEvAll
addpath('\\zserver\Code\Rigging\main'); %dat.expPath
addpath('\\zouni\Users\Experiment\Documents\MATLAB\CWanalysis\utility');
% addpath(genpath('C:\Users\Experiment\widefield_ns'));
addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\widefield_ns'));

addpath(genpath('C:\Users\Experiment\kilotrodeRig'));
addpath 'C:\Users\Experiment\npy-matlab'; %12/12/17

addpath 'C:\Users\Experiment\Documents\MATLAB\widefield\generalUtils';%LoadRotateInfo_2exps_SVD
addpath '\\basket\home\code\zouni20170127\wheelAnalysis';

run('C:\Users\Experiment\Documents\MATLAB\CWanalysis\utility\gndPaths.m');

expt.subject
expt.expDate
expt.expNum

%% load imaging data

nSV = 500;
figDestination = 'E:\gcamp_ns';

%trigInfo = loadTrigInfo(15); %18/1/18 ??
vidName_trig = '_corr_dFF'; %string 17/1/18

resizeFac2 = 0.25;
%     [expid, expid_passive, expid_base, ~, contrastList] ...
%         = expListCWpaper(expt.subject);

%     expList = [expid expid_passive];
%     for eee = 1:length(expList)
%         try
%expt = readExpsDatabase('ExpsDatabase_Monash.m', expList(eee));
thisDate = [expt.expDate(1:4) '-' expt.expDate(5:6) '-' expt.expDate(7:8)];
%expPath = dat.expPath(Exp.expt.subject, thisDate, Exp.iexp, 'main', 'master');
expPath = fullfile('\\zserver.cortexlab.net\Data\Subjects',...
    expt.expt.subject,thisDate,num2str(expt.iexp));

saveVpath = fullfile(figDestination, expt.subject, thisDate, num2str(expt.expNum));
mkdir(saveVpath);

%             try
%                 prepareAuxVids(Exp.expt.subject, thisDate, Exp.iexp);
%                 delete(fullfile(saveVpath, 'prepareAuxVids.png'));
%             catch err
%                 ng(eee) = err;
%                 disp(['ng prepareAuxVids:' num2str(expList(eee))]);
%             end


quickHemoCorrect(expPath, saveVpath, nSV);
%             newfig = mergefigs([1 2 3]);
%             set(newfig,'position',[1 1 1920 1200]);
%             screen2png(fullfile(saveVpath,'sanityCheckFigs'));

%             [U, V, t, mimg] = quickLoadUVt(expPath, nSV, saveVpath);
%              [newU, dV] = dffFromSVD(U, V, mimg);
%
%
%             writeUVtoNPY(newU, dV, fullfile(fileparts(saveVpath), 'svdSpatialComponents_corr_dFF'), ...
%                 fullfile(saveVpath, 'svdTemporalComponents_corr_dFF'));
%             copyfile(fullfile(saveVpath, 'svdTemporalComponents_corr.timestamps.npy'),...
%                 fullfile(saveVpath, 'svdTemporalComponents_corr_dFF.timestamps.npy'))
%             clear U mimg
%             close all


%                 %% purple
%                 %                 params.movieSuffix = 'purple';
%                 %                 [U, V, t, mimg] = quickLoadUVt(expPath, nSV, [], params);
%                 expRoot = fileparts(expPath);
%
%                 movieSuffix = 'purple';
%                 disp(['loading svdSpatialComponents_' movieSuffix '.npy'])
%                 U = readUfromNPY(fullfile(expRoot, ['svdSpatialComponents_' movieSuffix '.npy']), nSV);
%                 mimg = readNPY(fullfile(expRoot, ['meanImage_' movieSuffix '.npy']));
%                 disp(['loading svdTemporalComponents_' movieSuffix '.npy']);
%                 V = readVfromNPY(fullfile(expPath, ['svdTemporalComponents_' movieSuffix '.npy']), nSV);
%                 t = readNPY(fullfile(expPath, ['svdTemporalComponents_' movieSuffix '.timestamps.npy']));
%                 Fs = 1/mean(diff(t));
%
%                 %V = detrendAndFilt(V, Fs); 5/1/18 commented out
%
%                 if length(t)==size(V,2)+1 % happens if there was an extra blue frame at the end
%                     t = t(1:end-1);
%                 end
%                 [newU, dV] = dffFromSVD(U, V, mimg);
%
%                 writeUVtoNPY(newU, dV, fullfile(fileparts(saveVpath), 'svdSpatialComponents_purple_dFF'), ...
%                     fullfile(saveVpath, 'svdTemporalComponents_purple_dFF'));
%                 copyfile(fullfile(expPath, 'svdTemporalComponents_purple.timestamps.npy'),...
%                     fullfile(saveVpath, 'svdTemporalComponents_purple_dFF.timestamps.npy'))
%                 clear U V t mimg

%             end


%% load rotation info
%             if expid_base ~= expList(eee)
%                 Exp_base = readExpsDatabase('ExpsDatabase_SVDtest.m', expid_base);
%                 thisDate_base = [Exp_base.iseries(1:4) '-' Exp_base.iseries(5:6) '-' Exp_base.iseries(7:8)];
%                 expPath_base = dat.expPath(Exp_base.expt.subject, thisDate_base, Exp_base.iexp, 'main', 'master');
%                 Exps(1) = Exp; Exps(2) = Exp_base;
%                 rotateInfo = LoadRotateInfo_2exps_SVD({fileparts(expPath) fileparts(expPath_base)}, Exps);
%             end

%% correlation map
%pixelCorrelationViewerSVD(U, V)

%% a movie of the session
%quickMovieWithVids(mouseName, thisDate, expNum)

%% load timeline and block
%load(dat.expFilePath(Exp.expt.subject, thisDate, Exp.iexp, 'Timeline', 'master'));
load(dat.expFilePath(expt.expt.subject, thisDate, expt.expNum, 'block', 'local')); %10/6/18


%[~, ~, ~, ~, contrastList] = expListCWpaper(Exp.expt.subject); %4/8/17
contrast_cache = [];
for itr = 1:block.numCompletedTrials
    contrast_cache = [contrast_cache block.trial(itr).condition.visCueContrast];
end
contrastList = unique(contrast_cache);
[condInfo, cue, resp, stimContrast] = condInfoFromContrastList(contrastList);%10/8/17

%% select trials
trigInfo_c = trigInfo;
trigInfo_c.eventWindow = [-0.5 0.15]; %18/1/18
[okTrials, SVDNGTrials, saccadeTrials, blinkTrials, wheelMovTrials]...
    = oktrials(expt, trigInfo_c, 4*360, vidName_trig, vidName_trig);

TGTTR = [];
for icond = 1:length(condInfo)
    tgtTr = selecttrials(block.trial, cue(icond), resp(icond), ...
        trigInfo.resptime, trigInfo.nrep, stimContrast(icond));
    TGTTR{icond} = intersect(tgtTr, okTrials);
end

%             trigInfo = loadTrigInfo(11);
%             %[Savg_nf, trigPath] = trig_bhv_SVD_gcamp(Exp, trigInfo, vidName_trig, TGTTR, nSV);
%             [Savg_nf, trigPath] = trig_bhv_RealDerRec_gcamp(Exp, trigInfo, vidName_trig, TGTTR, resizeFac2);
%             movefile(trigPath, [trigPath '_DerRec_nfilt']);

%% rename/delete previous results
thisDate = ['20' num2str(expt.expDate(3:4)) '-' num2str(expt.expDate(5:6)) ...
    '-' num2str(expt.expDate(7:8))];
thisPath = fullfile(figDestination, expt.expt.subject, thisDate);
Vpath = [thisPath '\' num2str(expt.expNum) '\'];
savePath = fullfile(Vpath,['trig_' trigInfo.evoke '_' trigInfo.nrep ...
    '_' num2str(1000*trigInfo.eventWindow(1)) ...
    '-' num2str(1000*trigInfo.eventWindow(2))] );
if exist(savePath, 'dir')
    delete(fullfile(savePath,'*'));
    rmdir(savePath);
end

savePath = [savePath '_nfilt'];
if exist(savePath, 'dir')
    delete(fullfile(savePath,'*'));
    rmdir(savePath);
end

%[Savg_nf, trigPath] = trig_bhv_SVD_gcamp(Exp, trigInfo, vidName_trig, TGTTR, nSV);
[Savg, trigPath] = trig_bhv_RealDerRec_gcamp(expt, trigInfo, vidName_trig, TGTTR, resizeFac2);

%             [Savg, trigPath] = trig_bhv_SVD_gcamp(Exp, trigInfo, '_purple_dFF', TGTTR, nSV);
%             movefile(trigPath, [trigPath '_purple']);



%             %% sanity check
%             load(dat.expFilePath(Exp.expt.subject, thisDate, Exp.iexp, 'Timeline', 'master'));
%
%             tr = block.trial; tr = tr(1:block.numCompletedTrials);
%             cond = [tr.condition];
%             stimOn = [tr.stimulusCueStartedTime];
%             vcc = [cond.visCueContrast];
%             contrastLeft = vcc(1,:);
%             contrastRight = vcc(2,:);
%             choice = [tr.responseMadeID];
%             goCue = [tr.interactiveStartedTime];
%             responseTime = [tr.responseMadeTime];
%             reactionTime = [tr.responseMadeTime]-goCue;
%             feedback = [tr.feedbackType];
%             repNum = [cond.repeatNum];
%             trialStarts = [tr.trialStartedTime];
%             trialEnds = [tr.trialEndedTime];
%
%             [contrastCondDefs, ia, contrastCondsRaw] = unique(vcc', 'rows');
%
%             sw = block.stimWindowUpdateTimes;
%
%             %% get alignment between block and timeline
%             pd = Timeline.rawDAQData(:,2);
%             tt = Timeline.rawDAQTimestamps;
%
%             pdFlips = schmittTimes(tt,pd, [3 5]);
%
%             figure;
%             plot(sw, ones(size(sw)), '.');
%             hold on;
%             plot(pdFlips, ones(size(pdFlips))+1, '.');
%             ylim([-3 6])
%
%             blockToTL = makeCorrection(pdFlips(2:end-1), sw, false);
%
%             %% look at an event-triggered activity
%             % stim onset, just for single stimulus trials
%             inclTrials_ns = contrastRight==0 & contrastLeft==0;
%             %inclTrials_ns = contrastRight==0 | contrastLeft==0;
%
%             % load('E:\gcamp_ns\Theiler\2017-06-07\2\trig_stimCue_nrep_-500-800ms\infoVEv_Theiler_20170607_l_x_50.mat',...
%             %     'block_trial');
%             load('E:\gcamp_ns\Chain\2016-12-08\1\trig_stimCue_nrep_-500-800ms\infoVEv_Chain_20161208_o_x_0.mat',...
%                 'block_trial');
%             inclTrials_ds = logical(zeros(1, length(inclTrials_ns)));
%             inclTrials_ds([block_trial.trialNumber]) = 1;
%             %check 1
%             ismember(find(inclTrials_ds), find(inclTrials_ns))
%
%             eventTimes_ns = blockToTL(stimOn(inclTrials_ds));
%             eventTimes_ds = cell2mat(evokeTime_gcamp(Exp, trigInfo.evoke, find(inclTrials_ds), Timeline));
%             %check 2
%             isequal(eventTimes_ns, eventTimes_ds)
%
%             inclTrials = contrastRight==0 | contrastLeft==0;
%             eventTimes = blockToTL(stimOn(inclTrials));
%             eventLabels = contrastRight(inclTrials)-contrastLeft(inclTrials);
%             pixelTuningCurveViewerSVD(newU, dV, t, eventTimes, eventLabels, [-0.5 2.5])
%
%         catch err
%             ng_pipeline{eee,aaa} = err;
%             disp(err)
%         end
%
%
%     end

t = toc
%% Description

%{

Extract time series features and append to "(species)_(subject)_(channel)"
Apply initialization and save as "(species)_(subject)_(channel)_hctsa"

Run this after preprocess_kirill.m

%}

%% Settings
addDirPrefs_COS;
%preprocessSuffix = '_subtractMean_removeLineNoise';
dirPref = getpref('cosProject','dirPref');
%species = 'human';
%subject = '376';
%load_dir = fullfile(dirPref.rootDir, 'preprocessed',species,subject);

saveDir = '/home/daisuke/Documents/git/analysisImaging/wfAnalysis/hctsa_tmp';
loadName_combined = fullfile(saveDir, 'data_combined.mat');
hctsaName = fullfile(saveDir, 'hctsa.mat');



%% load channels to process
% load(fullfile(load_dir,['detectChannels_' subject]) ,'tgtChannels');

expt.subject = 'Confucious';
expt.expDate = '2026-04-25_1';
params.movieSuffix = 'amber';
params.useCorrected = 0;

data_proc = [];
for icond = 1:2
    switch icond
        case 1
            expt.expNum = 6;
            state{icond} = 'awake';
        case 2
            expt.expNum = 5;
            state{icond} = 'unconscious';
    end
    %savedata_prefix = sprintf('%s_%s_ch%03d', species, subject, thisCh);
    %loadName = fullfile(load_dir, [savedata_prefix preprocessSuffix '.mat']);
    p = ProtocolLoad_wf(expt.subject,expt.expDate,expt.expNum); %3/6/20

    thisDate = expt.expDate(1:10);
    thisSeries = str2num(expt.expDate(12:end));

    loadName = [dat.constructExpRef(expt.subject, thisDate, thisSeries, expt.expNum) '_'...
        p.xfile(1:end-2) '_' params.movieSuffix '_' num2str(params.useCorrected) '_singleTraces'];
    loaded = load(fullfile(saveDir,loadName), 'dur','mrec');
    data_proc(:,:,icond) = loaded.mrec;
end

%% preprocess_kirill.m
data = [];
data.data_raw = data_proc;
data.data_proc = data_proc;
data.preprocess_params = [];%params;
data.preprocess_string = [];%preprocess_string;
data.channel = 1;%thisCh;
%data.lobe = thisLobe;

%out_file = [savedata_prefix preprocess_string];
save(loadName_combined,'data', '-v7.3', '-nocompression');%'epochID'
clear data

%% main_hctsa_1_init.m Setup for HCTSA - training set
% Reformat into series x time matrix

% Training dataset
data_set = data_proc; % time x trials x conditions
% Get labels for each time-series
%   dimensions - (channels x trials x flies x conditions)
dims = size(data_set);
ids = cell(dims(2:end)); % details of each time-series
for tr = 1 : dims(2)
    for c = 1 : dims(3)
        ids{tr, c} = ['epoch:' num2str(tr) ',state:' state{c}];
    end
end

% Reformat to (series x time)
data_set = permute(data_set, [2 3 1]); % trials x conditions x time
data_set = reshape(data_set, [prod(dims(2:end)) dims(1)]); % Collapse all dimensions other than time
ids = reshape(ids, [prod(dims(2:end)) 1]); % Collapse labels also
% Create hctsa matrix
timeSeriesData = data_set;
labels = ids; % keywords are already unique
keywords = ids;
save(loadName_combined, 'timeSeriesData', 'labels', 'keywords','-append');
tic;
TS_Init(loadName_combined, 'hctsa', [false, false, false], hctsaName);
toc
disp('training set done');


%% main_hctsa_2.m
TS_Compute(true, [], [], [], hctsaName);


%% main_hctsa_3.m
hctsa = matfile(hctsaName, 'Writable', true);
TS_DataMat = hctsa.TS_DataMat;
TS_Quality = hctsa.TS_Quality;

% "Fatal" errors - treat as NaN
TS_DataMat(TS_Quality == 1) = NaN;
% Special value NaN
TS_DataMat(TS_Quality == 2) = NaN;
% Special value Inf
TS_DataMat(TS_Quality == 3) = Inf;
% Special value -Inf
TS_DataMat(TS_Quality == 4) = -Inf;
% Special value complex
TS_DataMat(TS_Quality == 5) = NaN;
% Special value empty
TS_DataMat(TS_Quality == 6) = NaN;


% % Check for other cases
% if any(TS_Quality(:) > 4)
%     tmp = unique(TS_Quality(:));
%     disp([file_string ' TS_Quality ' num2str(tmp)]);
% end

hctsa.TS_DataMat = TS_DataMat;
hctsa.TS_Quality = TS_Quality;


%%  below from main_hctsa_matrix.m
hctsa.valid_features = getValidFeatures(hctsa.TS_DataMat);

TS_Normalised = BF_NormalizeMatrix(hctsa.TS_DataMat, 'mixedSigmoid');
hctsa.TS_Normalised = TS_Normalised;

%% NMclassification_selectCh.m
load('/mnt/hctsa-market/COSproject/Stage 1/results_subtractMean_removeLineNoise/TS_Normalised_train_macaque_George_ch055_validate_human_376_ch134_accuracy.mat',...
    'order_f'); %order used in Stage 1

data_all = hctsa.TS_Normalised;%[trainData.(htcsaType); validateData.(htcsaType) ];
TimeSeries_all = hctsa.TimeSeries;%[trainData.TimeSeries; validateData.TimeSeries];
order_e{1}=1:20;%
order_e{2}=1:20;%
order_f_valid = intersect(order_f', find(hctsa.valid_features)','rows','stable');

% fig = showHCTSAbarcodes(data_all, TimeSeries_all, order_f, order_e);
condNames = {'awake','unconscious'};
nConds = numel(condNames);
condTrials = getCondTrials(TimeSeries_all, condNames);
data = [];
for icond = 1:nConds
    data{icond} = TS_Normalised(condTrials==icond, order_f);%_valid);
end
for icond = 1:nConds
    ax(icond)=subplot(2,1,icond);
    imagesc(data{icond}(order_e{icond},:));
    title(condNames{icond});
end
colormap(inferno);
linkcaxes(ax(:), [0 1]);
mcolorbar;
savePaperFigure(gcf,hctsaName(1:end-4));
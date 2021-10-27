function S = LoadOneStimulus(Cam, protocol, iStim, iTr, ResizeFactor, defineIndividual)
% Loads the desired repeats for one stimulus from raw data
%
% S = StackSet.LoadOneStimulus( Cam, p, iStim, iTr ) returns a stack with all
% the repeats for stimulus iStim, repeat number iTr. 
% Cam is an Imager object (in which Cam.DataDir specifies directory of the raw data)
% p is the  p-file protocol. It returns [] if the data can't be loaded.. iTr must be
% scalar
%
% S = StackSet.LoadOneStimulus( Cam, p, iStim, iTr, ResizeFactor )
% lets you specify the resize factor ResizeFactor <= 1 (DEFAULT: 1).
%
% S = StackSet.LoadOneStimulus( Cam, p, iStim, iTr, ResizeFactor, defineIndividual )
% if defineIndividual = true, find FrameRate, nFrames and Dur for each data
% (default:false)
%
% See also: tools.Imager, StackSet.LoadStacks

% 2013-10 Matteo Carandini
% 2014.3.3 DS added doIndividual option
% 2014.3.27 DS S.Values is in single precision

%TO DO: enable iTr to be a vector!!

persistent Data

if nargin<6
    defineIndividual = false;
end

if nargin<5 || isempty(ResizeFactor)
    ResizeFactor = 1;
end

if length(iTr) ~= 1
    error('iTr must be a scalar');
end

S = StackSet; % create an empty object

%%


DataDir = fullfile(...
    Cam.DataDir,[protocol.animal Cam.FileString],...
    num2str(protocol.iseries), num2str(protocol.iexp) );

if ~exist(DataDir,'dir'),
    fprintf('There is no directory %s. Exiting.\n',DataDir);
    S = [];
    return
end

%% preparation for reading the data files

if (isempty(Data) || ~strcmp(DataDir, Data.Dir)) && ~defineIndividual
    
    Data.Dir = DataDir;
    
    % sort the data files in ascending order
    
    DirContents = dir(sprintf('%s/*.mat',Data.Dir));
    Data.nFiles = length(DirContents);
    
    date = zeros(Data.nFiles,1);
    for ifile = 1:Data.nFiles
        date(ifile) = datenum(DirContents(ifile).date);
    end
    [~, Order] = sort(date);
    Data.FileNames = DirContents(Order);
    
    % figure out Data.nx, Data.ny, Data.nFrames, Data.Dur, Data.FrameRate, and Data.AllTimeStamps
    
    Data.AllTimeStamps = cell(Data.nFiles,1);
    NumTimeStamps = zeros(Data.nFiles,1);
    for iFile = 1:Data.nFiles
        fname = fullfile(Data.Dir, Data.FileNames(iFile).name);
        switch Cam.Type
            case 'PCO'
                [Data.nx,Data.ny,TimeStamps] = tools.LoadPCO(fname, true, false); %DS 3/3/14
            case 'PhotonFocus'
                [Data.nx,Data.ny,TimeStamps] = tools.LoadPhotonFocus(fname);
            otherwise
                error('We cannot read this camera type yet');
        end
        Data.AllTimeStamps{iFile} = TimeStamps;
        NumTimeStamps(iFile) = nnz(isfinite(Data.AllTimeStamps{iFile}));
    end
    
    Data.FrameRate = 1./nanmedian(diff(Data.AllTimeStamps{1}));
    Data.nFrames = floor(median(NumTimeStamps));
    Data.Dur = (Data.nFrames-1)/Data.FrameRate;

elseif defineIndividual
    Data.Dir = DataDir;
    
    DirContents = dir(sprintf('%s/*.mat',Data.Dir));
    Data.nFiles = length(DirContents);
    
    date = zeros(Data.nFiles,1);
    for ifile = 1:Data.nFiles
        date(ifile) = datenum(DirContents(ifile).date);
    end
    [~, Order] = sort(date);
    Data.FileNames = DirContents(Order);

    iseqnums = protocol.seqnums(iStim, iTr);
    fname = fullfile(Data.Dir, Data.FileNames(iseqnums).name);
    switch Cam.Type
        case 'PCO'
            [Data.nx,Data.ny,TimeStamps] = tools.LoadPCO(fname,true,false);
        case 'PhotonFocus'
            [Data.nx,Data.ny,TimeStamps] = tools.LoadPhotonFocus(fname);
        otherwise
            error('We cannot read this camera type yet');
    end
    Data.FrameRate = 1./nanmedian(diff(TimeStamps));
    Data.nFrames = nnz(isfinite(TimeStamps));
    Data.Dur = (Data.nFrames-1)/Data.FrameRate;
end


%% Assign

S.nFrames   = Data.nFrames;
S.FrameRate = Data.FrameRate;
S.TimeVec   = linspace(0,Data.Dur,S.nFrames); %unreliable
if defineIndividual
    S.TimeVec = TimeStamps(isfinite(TimeStamps));
end
S.ResizeFactor = ResizeFactor;
S.nRows = ceil(Data.nx * ResizeFactor);
S.nCols = ceil(Data.ny * ResizeFactor);

S.Info.animal = protocol.animal;
S.Info.iseries = protocol.iseries;
S.Info.iexp = protocol.iexp;
S.Info.nConds = 'Nrepeats'; % a string to differentiate between Nrepeats and Nstimuli for nCond;
S.Info.Stimulus = iStim;
S.Info.RepeatList = iTr;

S.nConds = 1;
S.PixelSize = Cam.MmPerPixel / ResizeFactor; %DS on 28.2.14
S.Description   = sprintf('%s-%d-%d-Stim %d iTr:%s',protocol.animal,...
    protocol.iseries,protocol.iexp,iStim, num2str(iTr));

S.Values = zeros(S.nRows, S.nCols, S.nFrames, S.nConds, 'single');%27/3/14 DS

%% reading the raw data files
iseqnums = protocol.seqnums(iStim, iTr);

fprintf('Loading stimulus %d, repeat %d\n',iStim,iTr);
% load stimulus iStim, repeat iTr

FileName = fullfile(Data.Dir, Data.FileNames(iseqnums).name);
switch Cam.Type
    case 'PCO'
        [~,~,TimeStamps,OneStack] = tools.LoadPCO(FileName, true, false); %DO NOT CORRECT MISSED FRAMES HERE
        S = AssignCondition(S, OneStack(:,:,isfinite(TimeStamps)), 1); % this resizes covert from uint32 to double
    case 'PhotonFocus'
        [~,~,TimeStamps,OneStack] = tools.LoadPhotonFocus(FileName);
        S = AssignCondition(S, OneStack(:,:,isfinite(TimeStamps)), 1); % this resizes
    case 'Dalsa'
        [~,~,TimeStamps,OneStack] = tools.LoadDalsa(FileName); % to be completed
    case '2photon'
        [~,~,TimeStamps,OneStack] = tools.Load2Photon(FileName); % to be completed
    otherwise
        disp('unknown camera type!');
end

S.nFrames = length(S.TimeVec);%is this ok for option 1 and 2?


fprintf('\n');

function [traces, t, mimgOriginalSize, roimap] = quickLoadS2Pt(expt, marketDir, protocol, params)
%[traces, t, mimg, roi] = quickLoadS2Pt(marketDir, protocol, params)
%
%OUTPUT
%traces{layers}[cellx,frames]: traces of ROIs or cells concatenated across repeats
%t{layers}[1,frames]: time in Thorsync time after reducing reduce thorsync frames to match thorimage frames
%mimgOriginalSize(x,y,layers): mean image after image registration
%roimap{layers}(x,y,cells):

%16/7/20 created

%TODO
%make this function faster. the bottole neck is the reading of thorsync
%file

modality = params.modality;
useCells = params.useCells;

%% detect 2p frame time in thorsync time
disp('grabbing 2p Frame Times');
[ frameTimes, TIdata ] = grabFrameTimes( expt, expt.expNum, 'all', marketDir); %[1 x 2pframes]
nlayers = size(frameTimes,2);
%frameTimes [frame x layer]

%frameTimes_rc{repeats}(frames, layer)
frameTimes_rc = parseFrameTimes(frameTimes, protocol.nrepeats, protocol.blockEndDelay);

%% load 2p traces
disp('load 2p traces')
nFramesTI = [];
traces = cell(nlayers,1);
roimap = cell(nlayers,1);
mimgOriginalSize = [];
for iplane = 1:nlayers
    s2pname = sprintf('%s/%s/%s/%d/F_%s_%s_plane%d_proc.mat', marketDir, ...
        expt.subject, expt.expDate, expt.expNum, expt.subject, expt.expDate, iplane);
    load(s2pname,'dat');
    
    if strcmp(modality, 'Fcell')
        traces_c = dat.Fcell{1};%[cells x 2p frames] %F
    elseif strcmp(modality, 'FcellNeu')
        traces_c = dat.FcellNeu{1};%[cells x 2p frames]
    elseif strcmp(modality, 'FcellCorrected')
        if ~isfield(dat,'Ff'); error('Ff not processed. Run add_deconvolution'); end
        traces_c = dat.Ff{1};
    elseif strcmp(modality, 'sp')
        if ~isfield(dat,'sp'); error('sp not processed. Run add_deconvolution'); end
        traces_c = dat.sp{1};
    end
    
    %% use only rois judged as cells
    if useCells
        cellIdx = find([dat.stat.iscell]); %only ROIs judged as cells
    else
        cellIdx = 1:length([dat.stat.iscell]); %all ROIs
    end
    traces{iplane} = [traces{iplane} traces_c(cellIdx,:)]; %[cells x frames]
    ncells = length(cellIdx);
    
    for irepeat = 1:protocol.nrepeats
        nFramesTI(iplane, irepeat) = dat.ops.fsroot{1}(irepeat).Nframes;
    end
    
    mimgOriginalSize(:,:,iplane) = mimgInFoV(dat);
    [roimap{iplane}] = ROI2D(dat, useCells);

end



%% reduce thorsync frames to match thorimage frames for each repeat
frameTimes_r = adjustFrameTimes_repeat(frameTimes_rc, nFramesTI); %frameTimes_r{repeat}(time,plane)


%% concatenate time across repeats ... t is different from frameTimes! 7/8/20
t = []; 
for irepeat = 1:protocol.nrepeats
    t = [t frameTimes_r{irepeat}']; %[planes x frames]
end


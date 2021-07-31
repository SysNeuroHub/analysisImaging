function infoInXml = loadThorlabsExperimentXml(xmlDirname, xmlFilename)
% loadThorlabsExperimentXml(xmlDirname, xmlFilename)
%
% 15/8/19 created from loadThorlabsExperiment. maybe not work for .raw
% format?

if nargin < 2
    xmlFilename = 'Experiment.xml';
end

verbose = 0;

imInfo = xml2struct(fullfile(xmlDirname, xmlFilename));

verString=imInfo.ThorImageExperiment.Software.Attributes.version;
verStringCell=strsplit(verString,'.');
verNum=cellfun(@str2num,verStringCell);

if verNum(1) < 3,
    if verbose
        warning(['This is a 2.4 or earlier version experiment. Version: ' verString ]);
    end
    verNumShort=2.4;
elseif verNum(1)==3 & verNum(2)==0,
    verNumShort=3;
    if verbose
        warning(['This is a 3.0 version experiment. Version: ' verString ]);
    end
elseif verNum(1)==3 & verNum(2)==1,
    verNumShort=3;
    if verbose
        warning(['This is a 3.1 version experiment. Version: ' verString ]);
    end
elseif verNum(1)==3 & verNum(2)==2,
    verNumShort=3;
    if verbose
        warning(['This is a 3.2 version experiment. Version: ' verString ]);
    end
else
    if verbose
        warning(['Tested to work with ThorImage 3.2.X and earlier data. May not work with this version: ' verString ]);
    end
    verNumShort=999;
end

%This part determines which channels were in use during the experiment for later file opening.
%There is extra code here that also shows how to create an array that give a simple boolean 'was this channel active' record.
%This boolean record is set up as: [ChanAActive? ChanBActive? ChanCActive? ChanDActive?]. So, if channels A and C are active, then 'activeChannels' will equal [1 0 1 0].
%This is not required for the rest of the function, but is an example of how to get that information for easier creation of other programs that might depend on which channels are active.
allChannelNames={'ChanA';'ChanB';'ChanC';'ChanD'};  %List of all channel names for use later
activeChannelNames={};  %List of only the channel names that were in use during the imaging experiment. Filled in below
activeChannels=zeros(4,1);  %Boolean record of active channels. Filled in below
numActiveChannels=length(imInfo.ThorImageExperiment.Wavelengths.Wavelength);  %Detect the number of active channels from the XML tag data.
zFastEnabled=str2num(imInfo.ThorImageExperiment.Streaming.Attributes.zFastEnable);
for I=1:numActiveChannels,    %Loop once for each actuve channel
    if numActiveChannels==1,
        activeChannelNames{I}=imInfo.ThorImageExperiment.Wavelengths.Wavelength.Attributes.name;  %Feed the string of the active channel name from deeply burried structure to an easier to use cell array
    else
        activeChannelNames{I}=imInfo.ThorImageExperiment.Wavelengths.Wavelength{I}.Attributes.name;  %Feed the string of the active channel name from deeply burried structure to an easier to use cell array
    end
    
    for J=1:4,  %Loop once for each channel that *might* be active
        if strcmp(allChannelNames{J},activeChannelNames{I}),  %Determine if the active channel name matches the channel name in slot J
            activeChannels(J) = 1;  %Set the boolean to 1 for easier scripting later
        end
    end
end

%Determine how many wells were scanned, according to the XML tags
if 0 & isfield(imInfo.ThorImageExperiment.Sample,'Wells'),
    numWells = ...
        str2num(imInfo.ThorImageExperiment.Sample.Wells.Attributes.rows) * ...
        str2num(imInfo.ThorImageExperiment.Sample.Wells.Attributes.columns);
    %Determine how many sub images were scanned, according to the XML tags.
    numSubImages = ...
        str2num(imInfo.ThorImageExperiment.Sample.Wells.SubImages.Attributes.subRows) * ...
        str2num(imInfo.ThorImageExperiment.Sample.Wells.SubImages.Attributes.subColumns);
else
    numWells=1;
    numSubImages=1;
end

%Determine how many Z slices were scanned, according to the XML tags.
%Version before 3.0 do not have ZEnable
if 0 == isfield(imInfo.ThorImageExperiment.ZStage.Attributes, 'enable'),
    zEnabled = 1;
else
    zEnabled = str2num(imInfo.ThorImageExperiment.ZStage.Attributes.enable);
end
zFastEnabled=str2num(imInfo.ThorImageExperiment.Streaming.Attributes.zFastEnable);
if 0 == str2num(imInfo.ThorImageExperiment.CaptureMode.Attributes.mode) & 1 == zEnabled,
    numZSlices = str2num(imInfo.ThorImageExperiment.ZStage.Attributes.steps);
elseif 1 == str2num(imInfo.ThorImageExperiment.CaptureMode.Attributes.mode) & 1 == zFastEnabled,
    numZSlices = str2num(imInfo.ThorImageExperiment.ZStage.Attributes.steps);
else
    numZSlices = 1;
end


%Determine how many time points were taken, according to the XML tags.
numTimepoints = str2num(imInfo.ThorImageExperiment.Timelapse.Attributes.timepoints);

%Determine how many X and Y pixels exist, according to the XML tags.
%Version before 3.0 do not have Modality
if 1 == isfield(imInfo.ThorImageExperiment,'Modality'),
    if 0 == str2num(imInfo.ThorImageExperiment.Modality.Attributes.primaryDetectorType),
        numXPix=str2num(imInfo.ThorImageExperiment.Camera.Attributes.width);
        numYPix=str2num(imInfo.ThorImageExperiment.Camera.Attributes.height);
        xyPixelSizeUM=imInfo.ThorImageExperiment.Camera.Attributes.pixelSizeUM;
        heightUM = num2str(times(numYPix, str2num(xyPixelSizeUM)));
        widthUM = num2str(times(numXPix, str2num(xyPixelSizeUM)));
    else
        numXPix=str2num(imInfo.ThorImageExperiment.LSM.Attributes.pixelX);
        numYPix=str2num(imInfo.ThorImageExperiment.LSM.Attributes.pixelY);
        xyPixelSizeUM=imInfo.ThorImageExperiment.LSM.Attributes.pixelSizeUM;
        heightUM=imInfo.ThorImageExperiment.LSM.Attributes.heightUM;
        widthUM=imInfo.ThorImageExperiment.LSM.Attributes.widthUM;
    end
else
    numXPix=str2num(imInfo.ThorImageExperiment.LSM.Attributes.pixelX);
    numYPix=str2num(imInfo.ThorImageExperiment.LSM.Attributes.pixelY);
    xyPixelSizeUM=imInfo.ThorImageExperiment.LSM.Attributes.pixelSizeUM;
    heightUM=imInfo.ThorImageExperiment.LSM.Attributes.heightUM;
    widthUM=imInfo.ThorImageExperiment.LSM.Attributes.widthUM;
end

%Determine a lot of other data about the image
numTimepoints=str2num(imInfo.ThorImageExperiment.Timelapse.Attributes.timepoints);
numFlybackFrames=str2num(imInfo.ThorImageExperiment.Streaming.Attributes.flybackFrames);
zStepSizeUM=imInfo.ThorImageExperiment.ZStage.Attributes.stepSizeUM;

if verbose
    disp([ 'Image size is ' num2str(numXPix) ' by ' num2str(numYPix) ' pixels.' ]);
    disp([ 'Data contains ' num2str(numZSlices) ' z slices.' ])
    disp([ 'Data has ' num2str(numTimepoints) ' time points.' ])
    disp([ 'Active channels: ' ])
    disp(activeChannelNames)
    disp([ 'Actual image size: width = ' widthUM ' um by height = ' heightUM ' um.' ])
    disp([ 'Pixel size is: ' xyPixelSizeUM ' um.' ])
    disp([ 'Z step size (if applicable) is: ' zStepSizeUM ' um.' ]);
end

infoInXml = [];
infoInXml.numXPix = numXPix;
infoInXml.numYPix = numYPix;
infoInXml.numZSlices = numZSlices;
infoInXml.numTimepoints = numTimepoints;
infoInXml.activeChannelNames = activeChannelNames;
infoInXml.widthUM = widthUM;
infoInXml.heightUM = heightUM;
infoInXml.xyPixelSizeUM = xyPixelSizeUM;
infoInXml.zStepSizeUM = zStepSizeUM;
infoInXml.numFlybackFrames = numFlybackFrames;
infoInXml.zEnabled = zEnabled;
infoInXml.numWells = numWells;
infoInXml.numSubImages = numSubImages;
infoInXml.zFastEnabled = zFastEnabled;
infoInXml.verNum = verNum;
infoInXml.verNumShort = verNumShort;

infoInXml.magnification = str2num(imInfo.ThorImageExperiment.Magnification.Attributes.mag);
infoInXml.frameRate = str2num(imInfo.ThorImageExperiment.LSM.Attributes.frameRate);

infoInXml.averageNum = str2num(imInfo.ThorImageExperiment.LSM.Attributes.averageNum); %2/4/20
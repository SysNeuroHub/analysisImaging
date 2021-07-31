function [FTOSignalValid, stimulusSignalValid] = filterFTObyStimulusFastZ(syncDataOut, infoInXml, showFigs)
% [FTOSignalValid, stimulusSignalValid] = filterFTObyStimulusFastZ(syncDataOut, infoInXml, showFigs)
% returns FTO of acquired frames, without flyback frames in fastz mode
% 10/2/20 DS created from filterFTObyStimulus
% 1/9/20 DS debugged captureRising not detected at recording onset

%All plotting commands are optional. They are in here to make it easier to grasp what is going on.
%IMPORTANT - Since the lines in Thorsync can be named anything, you will probably need to change a couple of lines to make the variable names consistent. Check the line right under the part called VARIABLE NAMING

%Introduction:
% In a perfect world, the number of frame trigger out rising edges would equal the number of recorded frames of data in a stimulus experiment. In reality, there are several places where this can vary.
% The too-long-didn't-read version is: this script determines which stimulus triggers are valid by looking at the capture active line.
% It determines which frame trigger outs are valid by looking at the validated stimulus triggers.

% Some more detail:
% 0. If the number of frames recorded in every stimulus is equal to the 'Max Frames' value set in Thorimage, then the frame trigger number should actual equal number of recorded frames. You do not need to use this script.
% 1. Since the stimulus trigger can fall at any time, the hardware does not know when to stop and take a little time to actually stop running the scanners. This will correspond to 1 - 2 extra frames scanned but not recorded. Frame triggers will be generated for these orphan frames, but there will be no data.
% 2. If the stimulus falls during a frame, the falling edge of that current frame will need to be kept since data for that current frame will be kept.
% 3. If the stimulus falls but then rises before the stimulus trigger (internal) has a chance to reset, then this is a false stimulus trigger and should be filtered out. This script does that.
% 4. Recording the capture active line is critical. Stimulus rising edges can only be detected while the capture active line is low. Therefore, this script determines an invalid stimulus trigger as one that rises while the capture active line is still high from a previousl stimulus trigger.

if nargin < 3
    showFigs = 0;
end

%syncDataOut=LoadSyncEpisodeFunction(fnam);

%VARIABLE NAMING - Because the naming in Thorsync can vary, set these lines here to correspond to the values where the 'Stimulus Trigger', 'Frame Trigger Out', and 'Capture Active' are actually being recorded.
Stimulus_Trigger_In=syncDataOut.frameTrigger;%Bleach_Complete;
Frame_Trigger_Out=syncDataOut.FrameOut;%Frame_Out;
Capture_Active=syncDataOut.CaptureActive;
time=syncDataOut.time;

%I don't actually do anything with imaging data in this script, but I left the lines in here commented as an example of how to load image data.
% fnam='c:\Users\cszymanski\Documents\data\stimulus trigger filtering\test_003\Experiment.xml'  %Substitute for the path to the imaging data.
% imdata=loadThorlabsExperimentRaw(fnam);  %In this case I left it raw when acquiring

if showFigs
    %Plot validated results
    figure(1);
    plot(time,Stimulus_Trigger_In,time,Frame_Trigger_Out,time,Capture_Active)
    title('Unvalidated signals')
    legend('Stimulus','FTO','Capture Active')
end

stimulusEvents=[0;int8(Stimulus_Trigger_In(2:end))-int8(Stimulus_Trigger_In(1:end-1))];  %Convert the stimulus signal to its numerical derivative. This way spikes will indicate events rather than rising or falling edges. Rising edges will be positive spikes while falling edges will be negative spikes. The '[0;' part is so that the array keeps the same length which makes things easier later. The 'int8(' parts are so that we do not end up with a logical data type. A logical data type cannot go negative so any falling edges would be lost.
if (Stimulus_Trigger_In(1)>0.5); stimulusEvents(1) = 1; end %7/2/20 happens when ThorSync is triggered by Stimulus_Trigger_In
FTOevents=[0;int8(Frame_Trigger_Out(2:end))-int8(Frame_Trigger_Out(1:end-1))];  %Take the numerical derivative of the frame trigger out signal.
captureEvents=[0;int8(Capture_Active(2:end))-int8(Capture_Active(1:end-1))];  %Take the numerical derivative of the capture active signal.


% plot(time,stimulusEvents,time,FTOevents)

%Create several arrays to hold the different event times.
FTORisingValid        =zeros(size(FTOevents));
FTOFallingValid       =zeros(size(FTOevents));
FTORisingInvalid      =zeros(size(FTOevents));
FTOFallingInvalid     =zeros(size(FTOevents));
stimulusRisingValid   =zeros(size(stimulusEvents));
stimulusFallingValid  =zeros(size(stimulusEvents));
stimulusRisingInvalid =zeros(size(stimulusEvents));
stimulusFallingInvalid=zeros(size(stimulusEvents));
captureRising         =zeros(size(captureEvents));
captureFalling        =zeros(size(captureEvents));




%All stimulus events are considered invalid until proven otherwise. See STIMULUS VALID/INVALID EXPLANATION below.
%We do however separate them into rising and falling edges.
stimulusRisingInvalid (find(stimulusEvents> 0.5))=1;
stimulusFallingInvalid(find(stimulusEvents<-0.5))=1;

%All capture active signal are automatically valid.
captureRising (find(captureEvents> 0.5))=1; %hey craig, this is not robust!
captureFalling(find(captureEvents<-0.5))=1; %hey craig, this is not robust!

captureRising(find((diff(captureRising) == 0) & (captureRising(1:end-1)==1))+1) = 0; %30/3/20
captureFalling(find((diff(captureFalling) == 0) & (captureFalling(1:end-1)==1))+1) = 0; %30/3/20

%if number of detected events are not even, add rising event at the beginning %1/9/20
if max(cumsum(captureRising)) < max(cumsum(captureFalling))
    captureRising(1) = 1;
end

%Plot some more
% figure(1); plot(time,Stimulus_Trigger_In,time,FTORisingValid  ,time,FTOFallingValid  )
% figure(2); plot(time,Stimulus_Trigger_In,time,FTORisingInvalid,time,FTOFallingInvalid)
% figure(3); plot(time,stimulusRising,time,stimulusFalling)

%There is the possibility of a stimulus signal arriving too soon after the previous trigger has fallen. In this case, the new stimulus rising and falling edge will be ignored. It is important to filter these out because of some unusual cases.
%STIMULUS VALID/INVALID EXPLANATION - The logic for this is that the first FALLING edge during a particular capture active high timeframe is valid. Also, the stimulus RISING edge immediately BEFORE it is valid. All other stimulus edges are invalid.

stimulusRisingInvalidIdx =find(stimulusRisingInvalid ==1);
stimulusFallingInvalidIdx=find(stimulusFallingInvalid==1);


%Implementing the stimulus explanation above.
for I=1:sum(captureRising),  %Do once for each capture active rising edge (which should equal the number of falling edges)
    captureRisingIdx =min(find((cumsum(captureRising ))==I));  %Do a cumsum on the capture rising array so that we can easily index to the I-th rising or falling event.
    captureFallingIdx=min(find((cumsum(captureFalling))==I));
    foundStimulusFallingIdx = min(find((stimulusFallingInvalidIdx > captureRisingIdx) ...
        & (stimulusFallingInvalidIdx < captureFallingIdx)));  %Find the first falling edge that is later than the rising edge and earlier than the falling edge.
    foundStimulusRisingIdx = max(find(stimulusRisingInvalidIdx < ...
        stimulusFallingInvalidIdx(foundStimulusFallingIdx)));  %Find the latest rising edge that is earlier than the falling that was just found.
    
    stimulusRisingValid   (stimulusRisingInvalidIdx (foundStimulusRisingIdx ))=1;  %Set the index for the found rising edge to high in the valid array
    stimulusFallingValid  (stimulusFallingInvalidIdx(foundStimulusFallingIdx))=1;  %Set the index for the found falling edge to high in the valid array
    stimulusRisingInvalid (stimulusRisingInvalidIdx (foundStimulusRisingIdx ))=0;  %Set the index for the found rising edge to low in the invalid array
    stimulusFallingInvalid(stimulusFallingInvalidIdx(foundStimulusFallingIdx))=0;  %Set the index for the found falling edge to low in the invalid array
end

stimulusFallingValid = [0; stimulusFallingValid(1:end-1)];%HACK 9/7/20
%Recreate a stimulus signal using only the properly detected stimulus signals
stimulusSignalValid=cumsum(stimulusRisingValid)-cumsum(stimulusFallingValid);

%The next section uses the newly validated stimulus triggers to filter out invalid frame triggers.

%For the first pass, find any FTO events that occur while the stimulus signal is high. These are automatically valid, real events.
%I use 0.5 as a threshold because these signals can be detected on either analog or digital channels. 
%Having a threshold of 0.5 for high vs. low signal will work for either of these cases.
FTORisingValid   (find(FTOevents> 0.5 & stimulusSignalValid>0.5))=1;
FTOFallingValid  (find(FTOevents<-0.5 & stimulusSignalValid>0.5))=1;
FTORisingInvalid (find(FTOevents> 0.5 & stimulusSignalValid<0.5))=1;
FTOFallingInvalid(find(FTOevents<-0.5 & stimulusSignalValid<0.5))=1;

%The indices of the events are easier to work with than the raw detected event stream for the next part. 
%Finding the events and putting their indices (time step number) in arrays for convenience.
%Dropping the valid/invalid part for the next part since we are only using the valid ones anyway
stimulusRisingIdx =find(stimulusRisingValid ==1);
stimulusFallingIdx=find(stimulusFallingValid==1);
FTOFallingInvalidIdx=find(FTOFallingInvalid==1);

%explanation:
%Because the stimulus may (and usually will) have a falling edge during a frame, 
%we need to include the FTO falling edge that comes immediately after the stimulus 
%falling edge in the final tally of valid FTO events.


%   ------------   <- Stimulus signal starts high
%              |   <- Stimulus falling edge occurring during a frame
%       -------+------------   <- FTO signal. Happens to be high when stimulus signal goes low
%       |      |           |
%       |<-  Frame Time  ->|   <- Falling edge after the stimulus falling edge.
%       |      |           |      Needs to be included since it is connected with a frame of data that was saved
%   ----|      ------------------


%Go through and determine which falling stimulus edges occur during frame acquisitions and which do not. 
%For the ones that fall during frames, find the nearest falling FTO edge after the stimulus falling edge and include it as 'valid'.
for I=1:sum(stimulusFallingValid),  %Iterate loop once for each stimulus falling edge that exists.
    if Frame_Trigger_Out(stimulusFallingIdx(I))>0.5,  %Check to see if that falling edge happens when the FTO data stream is high or now. If low, skip the following several lines. If high, proceed to look for the next FTO falling edge.
        tempStimulusFallingIdx=min(find((cumsum(stimulusFallingValid))==I));  %Find the index (time position) of the I-th stimulus falling edge.
        %Explanation: The cumsum command creates an array that is a running tally of the array. So after the first spike occurs, the cumsum value will be 1, after the second, it will be 2, etc.
        %We then find the first on in the cumcum array that has a value equal to I, which will be the exact time point of the falling edge.
        
        additionalFTOFallingEventIdx=min(find(FTOFallingInvalidIdx>tempStimulusFallingIdx));
        
        % additionalFTOFallingEvent=min(find((FTOFallingInvalidIdx-tempStimulusFallingIdx)>0));  %Take the difference of the falling edge time point and the event time points. The lowest resulting index that is still above zero will be the closest FTO falling edge event that is after the stimulus falling edge.
        FTOFallingValid  (FTOFallingInvalidIdx(additionalFTOFallingEventIdx))=1;  %Set the FTO edge valid to 1 at that time point to indicate it is valid.
        FTOFallingInvalid(FTOFallingInvalidIdx(additionalFTOFallingEventIdx))=0;  %Set the FTO edge invalid to 0 at that time point to indicate it is valid.
    end
    
    
% %         %make sure number of rising is equal to falling in each chunk
% %         FTORisingEvents = max(cumsum(FTORisingValid(1:stimulusFallingIdx(I))));
% %         FTOFallingEvents = max(cumsum(FTOFallingValid(1:stimulusFallingIdx(I))));
% %         if FTORisingEvents > FTOFallingEvents
% %             FTOFallingValid(stimulusFallingIdx(I))=1;
% %             FTOFallingInvalid(stimulusFallingIdx(I))=0;
% %         end
end

%Convert the event list back into a stream for the FTO. This valid one will not have the invalid, dropped signals.
FTOSignalValid=cumsum(FTORisingValid)-cumsum(FTOFallingValid);


%% fastz
if infoInXml.numZSlices > 1
    FTOeventsValid=[0;int8(FTOSignalValid(2:end))-int8(FTOSignalValid(1:end-1))];  %Take the numerical derivative of the frame trigger out signal.
    FTOSignalFinal = int8(zeros(size(FTOSignalValid)));
    for I=1:sum(captureRising),  %Do once for each capture active rising edge (which should equal the number of falling edges)
        
        thisCapture_Active = int8(zeros(size(Capture_Active)));
        captureRisingIdx =min(find((cumsum(captureRising ))==I));  %Do a cumsum on the capture rising array so that we can easily index to the I-th rising or falling event.
        captureFallingIdx=min(find((cumsum(captureFalling))==I));
        thisCapture_Active(captureRisingIdx:captureFallingIdx) = 1;
        
        thisFTOeventsValid = FTOeventsValid.*thisCapture_Active;
        
        frameID = 1:numel(find(thisFTOeventsValid > 0.5));
        frameID_woFB = [];
        movStride = infoInXml.numZSlices + infoInXml.numFlybackFrames;
        %ignore frames after maxFrame
        if (mod(numel(frameID), movStride) <= infoInXml.numZSlices) || (mod(numel(frameID), movStride) == 0)
            %if residual does not include FB, discard residual
            maxFrame = movStride * floor(numel(frameID)/movStride);
        else%if residual includes FB, keep the residual
            maxFrame = movStride * ceil(numel(frameID)/movStride);
        end
        for iplane = 1:infoInXml.numZSlices
            frameID_woFB = cat(2, frameID_woFB, iplane:movStride:maxFrame);
        end
        frameID_woFB = sort(frameID_woFB);
        frameID_discard = setxor(1:numel(frameID), frameID_woFB);
        
        FTORisingTime = find(thisFTOeventsValid > 0.5);
        thisFTOeventsValid(FTORisingTime(frameID_discard)) = 0;
        FTOFallingTime = find(thisFTOeventsValid < -0.5);
        thisFTOeventsValid(FTOFallingTime(frameID_discard)) = 0;
        FTOSignalFinal = FTOSignalFinal + cumsum(thisFTOeventsValid);
    end
    FTOSignalValid = FTOSignalFinal;
end

if showFigs
    %Plot validated results
    figure(2);
    plot(time,stimulusSignalValid,time,FTOSignalValid,time,Capture_Active);
    title('Validated signals only')
    legend('Stimulus','FTO','Capture Active')
end


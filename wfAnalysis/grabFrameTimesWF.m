function [ times ] = grabFrameTimesWF( expt, expr, frames)
% [ times ] = grabFrameTimesWF( expt, expr, frames)
% Gets frame times from widefield experiment
% Modified version of getFrameTimes.m by MK that doesn't take info struct
% TOBE FIXED: to use timeline

subject = expt.subject;
expDate = expt.expDate;

if nargin < 3
    frames = 'all';
end

timelineFile = strcat(expDate,'_',num2str(expr),'_',subject,'_','Timeline.mat');

try
    load(fullfile('\\zserver.cortexlab.net\Data\expInfo',subject,expDate,num2str(expr),timelineFile));
catch
    try
    load(fullfile('\\zserver.cortexlab.net\Data2\Subjects',subject,expDate,num2str(expr),timelineFile));
    catch
        load(fullfile('\\zubjects.cortexlab.net\Subjects',subject,expDate,num2str(expr),timelineFile));
    end
end


ind = strcmp({Timeline.hw.inputs.name}, 'neuralFrames');

nTotalFrames = Timeline.rawDAQData(end, ind);
TTLs = [0; diff(Timeline.rawDAQData(:, ind))];
idx = find(TTLs);

if isequal(frames, 'all')
    frames = 1:nTotalFrames;
else
    frames = frames(frames<=nTotalFrames);
end

times = Timeline.rawDAQTimestamps(idx(frames));

end




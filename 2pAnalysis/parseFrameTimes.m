function frameTimes_r = parseFrameTimes(frameTimes, nrepeats, dt_repeats)
%frameTimes_r = parseFrameTimes(frameTimes, nrepeats, dt_repeats)
%returns frame times in ThorSync of each repeat and layer {repeats}(frames, layer)
%frameTimes must be output of grabFrameTimes (frames x layer)
%14/2/20 created

if nargin < 3
    dt_repeats = 1; %gap between repeats [s]
end
if nargin < 2
    nrepeats = numel(find(diff(frameTimes(:,1)))) + 1;
end
nlayers = size(frameTimes,2);

% check number of repeats recorded in ThorSync %9/7/20
repeatEndIdx = find(diff(frameTimes(:,1)) > dt_repeats);
if nrepeats > length(repeatEndIdx)+1 
    warning(['parseFrameTimes: only ' num2str(length(repeatEndIdx)+1) 'repeats are recorded in ThorSync']);
    nrepeats = length(repeatEndIdx)+1;
end


%divide into repeats
frameTimes_r = cell(1,nrepeats);%9/7/20 was (nrepeats,1)
for ilayer = 1:nlayers
    repeatEndIdx = find(diff(frameTimes(:,ilayer)) > dt_repeats); 
    for irepeat = 1:nrepeats
        if irepeat == 1
            if nrepeats == 1
                theseIdx = 1:numel(frameTimes(:,ilayer));
            else
                theseIdx = 1:repeatEndIdx(1);
            end
        elseif (irepeat > 1) && (irepeat < nrepeats)
            theseIdx = repeatEndIdx(irepeat-1)+1:repeatEndIdx(irepeat);
        elseif irepeat == nrepeats
            theseIdx = repeatEndIdx(irepeat-1)+1:numel(frameTimes(:,ilayer));
        end
        frameTimes_r{irepeat}(:,ilayer) = frameTimes(theseIdx,ilayer);%[frame x layer]
    end
end
function frameTimes_r = adjustFrameTimes_repeat(frameTimes_rc, nFramesTI)
%frameTimes_r = adjustFrameTimes_repeat(frameTimes_rc, nFramesTI)
%frameTimes_r{repeat}(time,plane)

verbose = true;

nlayers = size(nFramesTI,1);
nrepeats = min(size(nFramesTI,2), size(frameTimes_rc,2));
for iplane = 1:nlayers
    for irepeat = 1:nrepeats
        if verbose
            disp(['L:' num2str(iplane) ' repeat:' num2str(irepeat) ...
                ' TI frames:' num2str(nFramesTI(iplane, irepeat)) ', TS frames:' num2str(length(frameTimes_rc{irepeat}))]);
        end
        nFramesTS = numel(frameTimes_rc{irepeat}(:,iplane));
        if nFramesTI(iplane, irepeat) ~= nFramesTS
            if nFramesTI(iplane, irepeat) == nFramesTS - 1
                frameTimes_r{irepeat}(:,iplane) = frameTimes_rc{irepeat}(1:end-1,iplane);
                if verbose
                    disp('TS frames adjusted according to heuristics');
                end
            else
                if verbose
                    disp('Unknown bug, cannot correct');
                end
            end
        else
             frameTimes_r{irepeat}(:,iplane) = frameTimes_rc{irepeat}(:,iplane); %do nothing
        end
    end
end
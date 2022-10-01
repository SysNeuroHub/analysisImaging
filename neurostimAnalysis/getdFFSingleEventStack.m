function singlePeriEventStack = getdFFSingleEventStack(singlePeriEventStack, baseWin, dFFmethod, doMedian)
% singlePeriEventStack = getdFFSingleEventStack(singlePeriEventStack, baseWin, dFFmethod, doMedian)


if ~doMedian
    switch dFFmethod
        case 1 %subtract by grand avg of prestimulus
            preMean = nanmean(nanmean(singlePeriEventStack(:,:,baseWin,:),3),4);
            singlePeriEventStack = 100*(singlePeriEventStack - preMean)./preMean; %dI/I [%]
            
        case 2 %subtract by each trial prestimulus
            preMean = nanmean(singlePeriEventStack(:,:,baseWin,:),3);
            singlePeriEventStack = 100*(singlePeriEventStack - preMean)./preMean; %dI/I [%]
            
        case 3 %subtract by each stimulus condition
            for icond = 1:nConds
                theseEvents = find(stimInfo.stimLabels == stimInfo.condLabels(icond));
                preMean = nanmean(nanmean(singlePeriEventStack(:,:,baseWin,theseEvents),3),4);
                singlePeriEventStack(:,:,:,theseEvents) = ...
                    100*(singlePeriEventStack(:,:,:,theseEvents) - preMean)./preMean; %dI/I [%]
            end
    end
else
    error('not yet implemented')
end

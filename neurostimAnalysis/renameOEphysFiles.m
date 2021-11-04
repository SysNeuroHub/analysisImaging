function [saveDir, saveFilename] = renameOEphysFiles(expt, DirBase_before)
%rename openEphys files
%
%[saveDir, saveFilename] = renameOEphysFiles(expt, DirFull_before)
% expt.subject = 'L4GCaMP6s_252';
% expt.expDate = '20210113';%'2021-01-13_1';
% expt.expNum = 4;
%
% DirBase_before: fullfile(rootDir,subject,'oephys','experiment1/recording15');

nExps = length(expt);

if ischar(DirBase_before)
    DirBase_before_tmp = DirBase_before;
    clear DirBase_before;
    DirBase_before{1} = DirBase_before_tmp;
end
% if ischar(DirBase_after)
%     DirBase_after_tmp = DirBase_after;
%     clear DirBase_after;
%     DirBase_after{1} = DirBase_after_tmp;
% end
for iexp = 1:nExps
    %     saveDir = fullfile(DirBase_after, expt(iexp).subject, 'processed',...
    %         expt(iexp).expDate, num2str(expt(iexp).expNum));
    saveDir = fileparts(expFilePathOEphys(expt(iexp).subject, expt(iexp).expDate,...
        num2str(expt(iexp).expNum)));
    saveFilename = ['oephysData_', expt(iexp).expDate, '_' num2str(expt(iexp).expNum)];
    
    if  exist(fullfile(saveDir, saveFilename) ,'dir')
        disp(['Already copied OEphys file ' DirBase_before{iexp} ' in ' fullfile(saveDir, saveFilename)]);
    else
        copyfile(DirBase_before{iexp}, fullfile(saveDir, saveFilename), 'f');
        disp(['Copied OEphys file ' DirBase_before{iexp} ' to ' fullfile(saveDir, saveFilename)]);
    end
end


%% for testing
% expt.subject = 'rat1gou';
% expt.expDate = '20211026';
% expt.expNum = 3;
%
% rootDir = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\MarmosetData';
%
% DirBase_before = fullfile(rootDir,expt.subject,'oephys','experiment1/recording15');
% DirBase_after = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\MarmosetData';
%
% renameOEphysFiles(expt, DirBase_before, DirBase_after)

function [fullpath, filename] = expFilePathOEphys(subject, expDate, expNum)
%[fullpath, filename] = expFilePathOEphys(subject, expDate, expNum)
%temporary function ... to be integrated with dat.expFilePath
%
%created from dat.expFilePath

DirBase = '\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\MarmosetData';
filename = 'structure.oebin';
fullpath = fullfile(DirBase, subject, expDate, 'processed', expDate, ...
    num2str(expNum), filename);
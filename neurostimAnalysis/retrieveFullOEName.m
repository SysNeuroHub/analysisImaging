function fullOEName = retrieveFullOEName(pth, partialOEName)
% Open ephys session folders are names as follows:
% 
%   [pth]/[prefix_]yyyy-mm-dd_HH-MM-SS[_suffix]/
% 
% On the marmolab rigs this means they should follow the naming convention:
%
%   [pth]/subject.paradigm.HHMMSS_yyyy-mm-dd_HH-MM-SS/
%
% Open Ephys GUI v0.5.5.? seems to have intruduced a bug whereby the open
% ephys data/time folder names are truncated, e.g.
%
%   [pth]/m1899.cuesaccade.HHMMSS_2022-03-11_HH-MM-SS/
%
% is truncated to:
%
%   [pth]/m1899.cuesaccade.HHMMSS/
%
% This function will scan the supplied path and the truncated directory name into an appropriate directory name.
%

currentDir = pwd;
assert(exist(pth,'dir') == 7,'%s is not a valid folder name or does not exist.',pth);

pat = ['^.*\', filesep, '(?<year>\d{4})\', filesep, '(?<month>\d{2})\', filesep, '(?<day>\d{2})\', filesep, '*$'];
pinfo = regexp(pth,pat,'names');

assert(~isempty(pinfo),'Supplied path should match: [blah]/yyyy/mm/dd/')

wd = cd(pth);

d = dir(fullfile(pth, [partialOEName '*']));

% keep only the folders/directories
d(~[d.isdir]) = [];

pat = '^(?<subject>.+)\.(?<paradigm>.+)\.(?<hour>\d{2})(?<minute>\d{2})(?<second>\d{2})$'; % <subject>.<paradigm>.<timestamp>
m = regexp({d.name},pat,'names');

ii = 1;
src = d(ii).name;
fullOEName = src;
% fullOEName = sprintf('%s_%s-%s-%s_%s-%s-%s', ...
%     d(ii).name,pinfo.year,pinfo.month,pinfo.day,m{ii}.hour,m{ii}.minute,m{ii}.second);

cd(currentDir);
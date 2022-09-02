% function dst = getFullOEName(pth)
% % dst = getFullOEName(pth)
% % created from fixoephys
% d.name = 'dum';
% 
% pat = ['^.*\', filesep, '(?<year>\d{4})\', filesep, '(?<month>\d{2})\', filesep, '(?<day>\d{2})\', filesep, '*$'];
% pinfo = regexp(pth,pat,'names');
% 
% pat = '^(?<subject>.+)\.(?<paradigm>.+)\.(?<hour>\d{2})(?<minute>\d{2})(?<second>\d{2})$'; % <subject>.<paradigm>.<timestamp>
% m = regexp({d.name},pat,'names');
% 
% ii=1;
% dst = sprintf('%s_%s-%s-%s_%s-%s-%s', ...
%     d.name,pinfo.year,pinfo.month,pinfo.day,m{ii}.hour,m{ii}.minute,m{ii}.second);


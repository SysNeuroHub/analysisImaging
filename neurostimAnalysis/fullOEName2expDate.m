function [expDate, m] = fullOEName2expDate(fullOEName)
%[expDate, m] = fullOEName2expDate(fullOEName)
pat = '^(?<subject>.+)\.(?<paradigm>.+)\.(?<hour>\d{2})(?<minute>\d{2})(?<second>\d{2})_(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})'; % <subject>.<paradigm>.<timestamp>
m = regexp(fullOEName,pat,'names');
expDate = sprintf('%s\\%s\\%s',m.year,m.month,m.day);

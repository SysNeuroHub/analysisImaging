function rawDir = getServerDirectory(animal)
% rawDir = getServerDirectory(animal)
% returns the directory of raw data 
%
% See also: getDirectory (returns directory of StackSet)

%2014-08-11 DS created
%2014-08-20 DS name changed from getRawDirectory

% TODO: animal name is not enough to discriminate when suffix is attached

% changed by SF on 2016-06-15 to update data directories following the
% merger of zserver and zserver2

directory = {'\\zraw.cortexlab.net\Data\VSFP\', '\\zraw.cortexlab.net\Data2\VSFP\', ...
    '\\zraw.cortexlab.net\Data3\VSFP\' '\\zserver.cortexlab.net\Data\GCAMP\' ...
    '\\zserver.cortexlab.net\Data\Intrinsic\' }; 

%OLD
% changed by MC on 2015-09-17 (to make it work from Rockefeller)
% directory = {'\\zraw.cortexlab.net\Data\VSFP\', '\\zraw.cortexlab.net\Data2\VSFP\', ...
%     '\\zraw.cortexlab.net\Data3\VSFP\' '\\zserver2.cortexlab.net\Data\GCAMP\' ...
%     '\\zserver2.cortexlab.net\Data\Intrinsic\' }; 

% OLD
% directory = {'\\zraw\Data\VSFP\', '\\zraw\Data2\VSFP\', '\\zraw\Data3\VSFP\' '\\zserver2\Data\GCAMP\' '\\zserver2\Data\Intrinsic\'}; 

for ii = 1:length(directory)
    dd = dir([directory{ii} animal '*']);
    if ~isempty(dd)
        rawDir = directory{ii};
        break
    end
end
if ~exist('rawDir', 'var')
    error(['Directory of raw data for ' animal ' was not found.']);
end

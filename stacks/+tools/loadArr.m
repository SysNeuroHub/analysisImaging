function [arr, meta] = loadArr(path)
% loadArr Fast load an array from a binary file
%   [arr, meta] = LOADARR(path) loads the numeric or character array saved
%   in <path>.bin and any associated meta data previously saved with it
%   (from <path>.mat). See also saveArr.
% 
%   Note that both files are required (<path>.bin and <path>.mat) to load
%   the array.
% 
% Part of Burgbox
%
% See also. tools.SaveArr, StackSet.LoadStacks, StackSet.LoadOneStimulus_ratio

% 2013-08 CB created

% load the info file
s = load([path '.mat']);

if nargout > 1
  if isfield(s, 'meta')
    meta = s.meta;
  else
    meta = [];
  end
end

%load the array data from binary file
binpath = [path '.bin'];
fid = fopen(binpath);
try
  arr = reshape(fread(fid, inf, ['*', s.arrPrecision]), s.arrSize);
  fclose(fid);
catch ex
  fclose(fid);
  rethrow(ex);
end

end
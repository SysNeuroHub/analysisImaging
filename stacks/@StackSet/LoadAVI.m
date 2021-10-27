function S = LoadAVI( FileName, RowFromTo, ColFromTo, FrameFromTo )

% LoadAVI loads the stacks from an AVI file.
%
% If the stacks have already been saved, it loads those. Otherwise, it
% loads them from the raw data.
%
%   S = LoadAVI( FileName )
%
%   S = LoadAVI( FileName, [RowFrom RowTo], [ColFrom ColTo], [FrameFrom FrameTo] )
%
% Example:
% S = StackSet.LoadAvi('\\ZSERVER\Code\MovieInspect\ExampleAVI.avi',[29 376],[113 471]);
% S.PlayCondition(1);
%
% See also: StackSet

% 2014-02 MC

% FileName = '\\ZSERVER\Code\MovieInspect\ExampleAVI.avi';

if nargin < 4, FrameFromTo = []; end
if nargin < 3,   ColFromTo = []; end
if nargin < 2,   RowFromTo = []; end

if isempty(  RowFromTo),   RowFromTo = [0 Inf]; end
if isempty(  ColFromTo),   ColFromTo = [0 Inf]; end
if isempty(FrameFromTo), FrameFromTo = [0 Inf]; end

if ~exist(FileName,'file')
    error('File does not exist');
end

oWaitBar = waitbar(0,'Loading...');

S = StackSet;

m = VideoReader(FileName);

S.FrameRate = m.FrameRate;
nFrames = m.FrameRate*m.Duration;

RowFrom  = max(        1,RowFromTo(1));
RowTo    = min(m.Height, RowFromTo(2));

ColFrom  = max(        1,ColFromTo(1));
ColTo    = min( m.Width, ColFromTo(2));

FrameFrom = max(        1, FrameFromTo(1));
FrameTo   = min(  nFrames, FrameFromTo(2));

S.nRows = RowTo-RowFrom+1;
S.nCols = ColTo-ColFrom+1;
S.nFrames = FrameTo-FrameFrom+1;

S.TimeVec = linspace(0,m.Duration,S.nFrames);

S.SpaceUnit = '';
S.nConds = 1;

S.Values = nan(S.nRows, S.nCols, S.nFrames, S.nConds);

for iFrame = 1:S.nFrames
    ThisFrame = mean(read(m,FrameFrom+iFrame-1),3);
    S.Values(:,:,iFrame,1) = ThisFrame( RowFrom:RowTo, ColFrom:ColTo );
    waitbar(iFrame/S.nFrames)
end

close(oWaitBar) 





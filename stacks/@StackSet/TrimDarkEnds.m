function S = TrimDarkEnds( S, mads, verbose )
% Trim in time to ensure mean across space are above threshold
%
% S = S.TrimDarkEnds(  )
% S = S.TrimDarkEnds( mads ) lets you specify how many MADs
% away to put the threshold (DEFAULT: 2)
% S = S.TrimDarkEnds( mads, verbose ) lets you select whether to show panel
% (default: true)
%
% See also: tools.fixMissingFrames (more complehensive tool to eliminate errorneous frames)

if nargin < 3
    verbose = true;
end

if nargin < 2
    mads = 2;
end

MeanSpaceAverage = mean(S.SpaceAverages,2);
Thresh = nanmedian( MeanSpaceAverage ) - mads* mad( MeanSpaceAverage(isfinite(MeanSpaceAverage)));

if verbose
    ThisFig = figure; clf;
    plot(S.TimeVec,MeanSpaceAverage(:),'k-'); hold on
    plot(S.TimeVec([1 end]),[Thresh Thresh],'r');
    set(gca,'xlim',[-0.1 S.Duration+0.1])
    xlabel('Time (s)');
    ylabel('Intensity');
    title('Should we trim?');
    
    answer = questdlg('Do you want to trim as indicated in the figure?','Trim?','Yes','No','No');
else
    answer = 'Yes';
end

switch answer
    case 'Yes'
        t1 = S.TimeVec( 1 + find(diff(MeanSpaceAverage > Thresh) ==  1) );
        t2 = S.TimeVec( 1 + find(diff(MeanSpaceAverage > Thresh) == -1) );
        if isempty(t1), t1 = S.TimeVec(1  ); end
        if isempty(t2), t2 = S.TimeVec(end); end
        S = S.Trim( [t1(1), t2(end)] ); % 2015-09-17 MC added (1) and (end)
end

if verbose
    close(ThisFig);
end

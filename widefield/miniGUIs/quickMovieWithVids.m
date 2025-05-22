function quickMovieWithVids(expt, params)
% quickMovieWithVids(mouseName, thisDate, expNum)

%reposName = 'master';
nSV = 250;

mouseName = expt.subject;
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(mouseName, thisDate, thisSeries, expt.expNum, 'widefield','master'));

%expPath = dat.expPath(mouseName, thisDate, expNum, 'main', 'master');

[U, V, t] = quickLoadUVt(expPath, nSV, [], params);

load(dat.expFilePath(mouseName, thisDate, thisSeries, expt.expNum, 'Timeline', 'master'));
traces = prepareTimelineTraces(Timeline);

auxVid = prepareAuxVids(mouseName, thisDate, expNum);

writeMovieLocation = fullfile(expPath, sprintf('widefield_%s_%s_%d', mouseName, thisDate, expNum));
movieWithTracesSVD(U, V, t, traces, writeMovieLocation, auxVid);

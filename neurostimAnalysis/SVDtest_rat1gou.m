expt.subject = 'rat1gou';
expt.expDate = '20211026';
expt.expNum = 2;


%% SVD
nSV = 200;%1000;
thisDate = expt.expDate;
thisSeries = 1;
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;

params.movieSuffix = 'green';% 'purple''corr_dFF'; %'blue'  %U for corr_dFF can have NANs..
params.useCorrected = 0;

disp('Loading widefield data');
disp(expt)
[U, V, t, mimg, mask] = quickLoadUVt(expPath, nSV, saveVpath, params);
Fs = 1/median(diff(t));

movieWithTracesSVD(U,V,t)
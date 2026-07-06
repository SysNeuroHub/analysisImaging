
load('/mnt/dshi0006_market/Subjects/Confucious/2026-02-27_3/dataSummary_amber.mat');
y = dataSummary.imageMeans;
[a,b,yfit] = fitExp1D(y);
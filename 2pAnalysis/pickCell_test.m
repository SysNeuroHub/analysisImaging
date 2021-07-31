

% load('M:\Subjects\232\2020-07-16_1\7\F_232_2020-07-16_1_plane1_proc.mat','dat');
load('\\storage.erc.monash.edu.au\shares\MNHS-dshi0006\Subjects\97138\2020-08-03_1\2\F_97138_2020-08-03_1_plane1_proc.mat')

xMapmn = xMapm / nx;
pickCell_test(dat, xMapmn); 
%TODO pickCell_test
%CANCELLED make redreq-figure and redraw mean image standalone?
%link the ichosen between the two (original new_main does it)
%DONE weird bug of image resizing first time clicking
%show colorbar
%DONE show cell property as text
%currently 2nd input must be [0 1]. How can it be original?
%DONT use dat instead mimg and iclust
%should be compatible all ROIs / only cells
 
% fields in h:
% mimg
% xlim
% ylim
% axis (renamed from axes2)
% stat.ipix (not 2d ROI map)
% cl.Ly
% cl.Lx
% cl.rands will be replaced by some response property of every ROI
% return cellID??


%% test without using pickCell_test
dat_selected = initializeAfterLoading_test(dat);
figure;
ax=gca;
[I, ichosen] = redraw_meanimg_test(ax,dat_selected);
%ichosen is not updated by clicking...


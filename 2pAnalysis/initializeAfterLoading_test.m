function dat_new = initializeAfterLoading_test(dat)
%h = initializeAfterLoading(h, loadname)
%apply audomatic ROI labeling & returns h suitable for manual ROI labeling in new_main.m
%
% h = initializeAfterLoading(h, loadname, skipClassify)
% if skipClassify = true, dont apply automatic ROI labeling
%
% if the user selected a file, do all the initializations


dat_new.stat = dat.stat(find([dat.stat.iscell])); %4/8/20 use only ROI judeged as cells
% if h.init==0
%     h.dat.filename = loadname;%fullfile(filepath1, filename1);
dat_new.cl.Ly       = numel(dat.ops.yrange);
dat_new.cl.Lx       = numel(dat.ops.xrange);

% make up iclut here
[dat_new.res.iclust, dat_new.res.lambda, dat_new.res.lambda0] =...
    getIclust_test(dat.stat, dat.cl);
dat_new.res.iclust = reshape(dat_new.res.iclust, dat_new.cl.Ly, dat_new.cl.Lx);
% %     h.dat.res.lambda = reshape(h.dat.res.lambda, h.dat.cl.Ly, h.dat.cl.Lx);
%

% some functional property of cells. will be replaced
dat.ops.Nk = numel(dat.stat);
rands_orig   = .1 + .65 * rand(1, dat.ops.Nk);
dat_new.cl.rands        = rands_orig;


dat_new.ylim = [0 dat_new.cl.Ly];
dat_new.xlim = [0 dat_new.cl.Lx];


dat_new.ichosen = 1;

ops = dat.ops;
if isfield(ops, 'mimg1') && ~isempty(ops.mimg1)
    dat_new.mimg = ops.mimg1(ops.yrange, ops.xrange);
    %     dat_new.mimg(:,:,dat_new.maxmap) = ops.mimg1(ops.yrange, ops.xrange);
    %     h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end

pathtoAllenBrain = '/home/daisuke/Dropbox/GSlab/MRI_alignment/Kim2023_test/pattern_generation/Brain_template.nii';
saveName = '/home/daisuke/Dropbox/GSlab/MRI_alignment/Kim2023_test/matlab_functions/Allen_pills_mask';

b = niftiread(pathtoAllenBrain); %allen standard brain 
bi =niftiinfo(pathtoAllenBrain);
%v=niftiread('/home/daisuke/Dropbox/GSlab/MRI_alignment/Kim2023_test/tmpD/T2w_resample_in_Allen.nii.gz');


p = b;
bi.Filename = '';
se = strel("sphere", 10); %13

%% pills left of the brain
p0_l=0*p;
p0_l(10:12,70,45:90)=1; %lateral, height, AP
p0_l = imdilate(p0_l, se);

%% pills right of the brain
p0_r = flipud(p0_l);

%% pills front of the brain
p0_f = 0*p;
p0_f(45:69,73,120)=1;
p0_f = imdilate(p0_f,se);

p = p0_f + 2*p0_r + 3*p0_l;

se2 = strel("sphere", 10);
bd = imdilate(b, se2);

p = p.*(~bd); %exclude brain region

niftiwrite(p, saveName, bi);
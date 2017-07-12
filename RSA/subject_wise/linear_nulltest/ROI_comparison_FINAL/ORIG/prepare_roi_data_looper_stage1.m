clear all;close all;clc;

warning('on','all')

data_path{1}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session1/post_processed_normal/bramila/';
data_path{2}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session2/post_processed_normal/bramila/';
data_path{3}='/triton/becs/scratch/braindata/kauttoj2/Memento/2015/preprocessed/Session3/post_processed_normal/bramila/';

HOME = pwd;

S = {'S7','S8', 'S9', 'S10','S12','S13','S15','S16','S17','S19','S21','S22','S23'}; % S13 incomplete

%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%S = S(1:5);
%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MAX_IMAGES = [1369,1369,1369];
MOVIE_LENGTHS = [35*60+20,35*60+17,35*60+8];
N_max_vols = 700;
TR = 1.56;
nworker = 5;
searchlightRad_mm = 7;
voxSize_mm=[2,2,2];

%%-----------------------------------------

toolboxRoot = '/triton/becs/scratch/braindata/kauttoj2/code/RSAtoolbox';
addpath(genpath(toolboxRoot));
BRAMILAPATH = '/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila';
addpath(BRAMILAPATH);
addpath('/triton/becs/scratch/braindata/kauttoj2/code/connISC_toolbox');
addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/RSA');

nii=load_nii('../ORIG_KRON_mask.nii');
EPI_mask=nii.img;
fprintf('EPI mask has %i voxels\n',nnz(EPI_mask));

nii = load_nii('/triton/becs/scratch/braindata/kauttoj2/code/RSAtoolbox/Templates/grey.nii');
PATTERN_mask = EPI_mask.*(nii.img>0.10);

maskID=find(PATTERN_mask>0);
save_nii_oma(PATTERN_mask,'PATTERN_mask.nii');
fprintf('PATTERN mask has %i voxels\n',nnz(PATTERN_mask));

% nii=load_nii('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/RSA/subject_wise/loops/radius5/meanmap_kauppi_p005_extent50.nii');
% ROI_mask = nii.img.*PATTERN_mask;
% ROI_mask=ROI_mask>0;
% PRECUNEUS_mask = 0*ROI_mask;
% PRECUNEUS_mask(51+(-2:2),30+(-2:2),54+(-2:2))=1;
% PRECUNEUS_mask=(PRECUNEUS_mask.*PATTERN_mask)>0;
% ANGULAR_mask = 0*ROI_mask;
% ANGULAR_mask(67+(-2:2),35+(-2:2),61+(-2:2))=1;
% ANGULAR_mask=(ANGULAR_mask.*PATTERN_mask)>0;
%my_masks = {ROI_mask,PRECUNEUS_mask,ANGULAR_mask};

% Other data
rad_vox=searchlightRad_mm./voxSize_mm;
minMargin_vox=floor(rad_vox);
% create spherical multivariate searchlight
[x,y,z]=meshgrid(-minMargin_vox(1):minMargin_vox(1),-minMargin_vox(2):minMargin_vox(2),-minMargin_vox(3):minMargin_vox(3));
sphere=((x*voxSize_mm(1)).^2+(y*voxSize_mm(2)).^2+(z*voxSize_mm(3)).^2)<=(searchlightRad_mm^2);  % volume with sphere voxels marked 1 and the outside 0
sphereSize_vox=[size(sphere),ones(1,3-ndims(sphere))]; % enforce 3D (matlab stupidly autosqueezes trailing singleton dimensions to 2D, try: ndims(ones(1,1,1)). )
% compute center-relative sphere SUBindices
[sphereSUBx,sphereSUBy,sphereSUBz]=ind2sub(sphereSize_vox,find(sphere)); % (SUB)indices pointing to sphere voxels
sphereSUBs=[sphereSUBx,sphereSUBy,sphereSUBz];
ctrSUB=sphereSize_vox/2+[.5 .5 .5]; % (c)en(t)e(r) position (sphere necessarily has odd number of voxels in each dimension)
ctrRelSphereSUBs=sphereSUBs-ones(size(sphereSUBs,1),1)*ctrSUB; % (c)en(t)e(r)-relative sphere-voxel (SUB)indices
nSearchlightVox=size(sphereSUBs,1);

roifile = [BRAMILAPATH,filesep,'external',filesep,'rois_Power264.mat'];
load(roifile);

selected_rois = false(1,length(rois));
all_mask = 0*PATTERN_mask;
k=0;
for i=1:length(rois)
    center = rois(i).centroid;
    if norm(ind2mni(center)-rois(i).centroidMNI)>2.01
        error('bad coordinate!')
    end    
    mask = 0*PATTERN_mask;
    for j=1:size(ctrRelSphereSUBs,1)
        mask(ctrRelSphereSUBs(j,1)+center(1),ctrRelSphereSUBs(j,2)+center(2),ctrRelSphereSUBs(j,3)+center(3))=1;
    end
    mask=mask.*PATTERN_mask;
    if nnz(mask)>=0.5*nSearchlightVox
        k=k+1;
        my_masks{k}=find(mask);        
        selected_rois(i)=true;
        all_mask(mask>0)=k;                        
    end
    clear mask;
end    
my_rois = rois(selected_rois);

cd(HOME);
save_nii_oma(all_mask,'all_masks_numbered.nii');
save('roi_and_mask_data','my_masks','my_rois','PATTERN_mask','all_mask','searchlightRad_mm','ctrRelSphereSUBs','sphereSUBs');

fprintf('Using total %i masks\n',length(my_masks));
%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%my_masks=my_masks(122);
%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
clear nii ROI_mask all_mask

for ses=1:3
    for subj_ind = 1:length(S)
        subj = S{subj_ind};                       
        cd(data_path{ses});        
        M = get_nii_frame([subj,'_mask_detrend_fullreg_filtered_zscored.nii']);
        if M<MAX_IMAGES(ses)
            MAX_IMAGES(ses)=M;
            disp(M)
            warning('ses %i: %s has less volumes than expected!',ses,subj);
        end
    end
end
cd(HOME);

% create local cluster
try
    myCluster = gcp('nocreate');
    if isempty(myCluster)
        %delete(gcp)
        myCluster = parcluster('local');
        myCluster.NumWorkers=nworker;
        parpool(myCluster);
    end
    N_workers = myCluster.NumWorkers;
catch err % old matlab?
    if ~matlabpool('size')
       eval(['matlabpool local ',num2str(nworker)]);
    end
    N_workers = matlabpool('size');
end

fprintf('\n-- loading and saving data ---\n')
parfor subj_ind = 1:length(S)
    subj = S{subj_ind};
    fprintf('...Subject %s\n',subj);        
    all_roi_data=cell(1,3);
    for ses=1:3        
        cd(data_path{ses});
        data = load_nii_mask([subj,'_mask_detrend_fullreg_filtered_zscored.nii'],my_masks,N_max_vols);        
        for k=1:length(data)
            if nnz(isnan(data{k}))>0
                error('NaN''s found!!!!')
            end
        end                
        all_roi_data{ses}=data;
    end    
    cd(HOME);
    res = save_data(all_roi_data,[subj,'_roi_data.mat']);
    if res == 0
        error('Failed to save!')            
    end
end
cd(HOME);

fprintf('\n ALL DONE!\n');




